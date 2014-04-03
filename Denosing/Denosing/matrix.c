#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>


/*这里有个内存泄露的问题*/
double* dotSub( double* iAry1, double* iAry2, int iRows, int iCols )
{
	int i,j;
	double* result;
	result = ( double* )malloc( iRows * iCols *sizeof( double ) );
	for( i = 0; i < iRows; i++)
		for( j = 0; j < iCols; j++)
			*( result + i * iCols + j ) = *( iAry1 + i * iCols + j ) - *( iAry2 + i * iCols + j );
	return result;
}


double* dotMulti( double* iAry1, double* iAry2, int iRows, int iCols )
{
	int i,j;
	double* result;
	result = ( double* )malloc( iRows * iCols *sizeof( double ) );
	for( i = 0; i < iRows; i++)
		for( j = 0; j < iCols; j++)
			*( result + i * iCols + j ) = *( iAry1 + i * iCols + j ) * (*( iAry2 + i * iCols + j ));
	return result;
}


double* dotDiv( double* iAry1, double* iAry2, int iRows, int iCols )
{
	int i,j;
	double* result;
	result = ( double* )malloc( iRows * iCols *sizeof( double ) );
	for( i = 0; i < iRows; i++)
		for( j = 0; j < iCols; j++)
			*( result + i * iCols + j ) = *( iAry1 + i * iCols + j ) / *( iAry2 + i * iCols + j );
	return result;
}


double mean( double* iAry,int iRows, int iCols )
{
	int i,j;
	double sum  =0.0;
	for( i = 0; i < iRows; i++ )
		for( j = 0; j <iCols; j++){
			if( *( iAry + i * iCols + j ) > 0 ){
				sum +=  *( iAry +i * iCols + j );
			}else{
				sum +=  -*( iAry +i * iCols + j );
			}
		}
	sum /= iRows * iCols;
	return sum;
}