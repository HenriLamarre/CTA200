#include <stdio.h>
#include "average.h"


int main() {
    double arr[] = {1.0, 2.0, 3.0, 4.0};

    double result = average(4, arr);
    int n;
    double var;
    var = 0;    	
    for(n=0;n<4;n++){
	var += (arr[n] - result) * (arr[n]-result);
	}
    var /= 4;

    printf("The average of 1, 2, 3 and 4 is: %.4f\n", result);
    printf("The variance is: %.4f\n", var);
    return 0;    
/* Added a computation for the variance of the array */
}

