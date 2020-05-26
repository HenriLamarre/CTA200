import numpy as np
import matplotlib.pyplot as pl
import trident
import yt
from yt.units import kpc,mp
from yt.fields.particle_fields import add_volume_weighted_smoothed_field

def _ElectronDensity(field, data):
    xh= 1. - data['PartType0','Metallicity_00'] -\
        data['PartType0','Metallicity_01']
    ne= xh * data['PartType0','Density'] *\
        data['PartType0','ElectronAbundance'] / mp
    return ne

def find_normal_vector(ang_mom_sphere):
    ''' 
    Finds the galaxy's normal vector by averaging the angular momentas
    of the particles
    '''
    angular_momentum =[np.mean(ang_mom_sphere['angular_momentum_x']),
                       np.mean(ang_mom_sphere['angular_momentum_y']),
                       np.mean(ang_mom_sphere['angular_momentum_z'])]
    normalized_ang = np.array(angular_momentum)/np.linalg.norm(angular_momentum)
    return normalized_and

def load_data(ray_start, ray_direction, low_res = True):
    if low_res:
        # Import some FIRE2 data
        fname = '/mnt/raid-project/murray/lakhlani/FIRE2_core/'+\
            'm12c_res56000/output/snapshot_600.hdf5'
        ds = yt.load(fname)
        # The center of the Halo was recovered in another notebook
        center_Halo = [25277.66673046,34505.21241664,32868.48520185]
    sp = ds.sphere(center_Halo, (20, "kpc"))

    # Put the earth about 8kpc from the 
    # center where the box units are in kpc/h
    

    z_axis = np.array([ 0.90088129,-0.07389156,-0.42772998])
    x_axis = yt.ortho_find(z_axis)[1]
    y_axis = yt.ortho_find(z_axis)[2]
    
    ray_start_conv = x_axis*ray_start[0]+y_axis*ray_start[1]+\
        z_axis*ray_start[2]
    ray_direction_conv = x_axis*ray_direction[0]+y_axis*ray_direction[1]+\
        z_axis*ray_direction[2]

    # Adds the field to the data set
    ds.add_field(('PartType0', '_ElectronDensity'), 
                 function = _ElectronDensity, 
                 sampling_type = 'particle', units='1/cm**3')

    # smooths a particle field into a continuous field
    fn = add_volume_weighted_smoothed_field("PartType0", 
                                            "particle_position",
                                            "particle_mass",
                                            "smoothing_length","density",
                                            "_ElectronDensity",
                                            ds.field_info)

    # Labels the new field
    ds.field_info.alias(('gas', '_ElectronDensity'), fn[0])
    ds.derived_field_list.append(('gas', '_ElectronDensity'))

    return ds, center_Halo, sp, ray_direction_conv, ray_start_conv

def make_light_ray_plots(filename, plotting = False,
                         ray_start=[0,0,0], ray_direction=[1,0,0]):
    ''' 
    This function computes the dispersion measure along a light ray
    going through a simulated galaxy.
    ---------Arguments------------
    plotting: If Galaxy plots should be made

    ray_start: The position of the start of the ray
    with respect to the center of the Halo

    ray_direction: The position of the end of the ray
    with respect to the start of the ray
    
    '''
    
    ds, center_Halo, sp, ray_direction_conv, ray_start_conv =\
        load_data(ray_start, ray_direction)
    
    # RAY CREATION
    ray_end = list(center_Halo + ray_start_conv + ray_direction_conv)
    ray = trident.make_simple_ray(ds,
                                 start_position=ray_start_conv+center_Halo,
                                 end_position=ray_end,
                                 fields = ['_ElectronDensity','density'],
                                 data_filename="{}.h5".format(filename))
    
    # GALAXY PLOTS
    if plotting:
        # Creates three on-axis projection plots
        axes=['x','y','z']
        for axis in axes:
            prj_x = yt.ProjectionPlot(ds, axis, 'density', 
                                      width=10*kpc, data_source = sp,
                                      center=center_Halo)
            prj_x.annotate_ray(ray)
            prj_x.save(filename+'_galaxy_{}.pdf'.format(str(axis)))
        
    # COMPUTING THE DM
    length = np.linalg.norm(np.array(ray_end)-np.array(center_Halo)-\
                                np.array(ray_start_conv))
    int_elec_dens = []
    int_total = 0
    # Finding the entry where arc_length>1kpc:
    for i in range(len(ray.data['_ElectronDensity'])):
        # Convert in pc/cm^-3 from kpc/h
        int_total += ray.data['_ElectronDensity'][i]*\
            ray.data['dl'][i]*60000*1000/0.7
        int_elec_dens.append(float(int_total))
    print(ray.data['dl'])
        
    # LIGHT RAY DATA PLOTS
    fig, axs = pl.subplots(3)
    
    axs[0].plot(np.linspace(0,length,len(ray.data['density'])), 
                np.log10(ray.data['density']))
    axs[0].set_xlabel('Arc Length (kpc/h)')
    axs[0].set_ylabel('Density '+r'$\log\left(\frac{g}{cm^3}\right)$')


    axs[1].plot(np.linspace(0,length,len(ray.data['_ElectronDensity'])), 
                np.log10(ray.data['_ElectronDensity']))
    axs[1].set_xlabel('Arc Length (kpc/h)')
    axs[1].set_ylabel('Electron Density '+r'$\log(cm^{-3})$')
        
    axs[2].title.set_text('DM = {}'.format(round(int_elec_dens[-1], 2)))
    axs[2].plot(np.linspace(0,length,len(int_elec_dens)), 
                      np.log10(int_elec_dens))
    axs[2].set_xlabel('Arc Length (kpc/h)')
    axs[2].set_ylabel('DM '+r'$\log(cm^{-2})$')

    pl.subplots_adjust(left = 0.20,right = 0.9,bottom = 0.2,top = 0.9,
                           wspace = 0.4,hspace = 1.5)
    fig.savefig(filename+'_subplots.pdf')
    return int_elec_dens[-1]
