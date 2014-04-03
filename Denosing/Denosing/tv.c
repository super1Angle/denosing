#include "tv.h"
#include <math.h>


int iter;   //迭代次数
double dt;   //时间步长
int ep;    //加入的常数值，防止奇异方差

//声明一些偏微分
/***************************************************/
//int flag = 0;
//int i, j;

void tv( Img* I, Img* I0, double lam)
{
	int i, j, m;
	double *tmp = ( double* )malloc( width * height  * sizeof( double ) );
	double *tmp1 = ( double* )malloc( width * height  * sizeof( double ) );

	double *dx = ( double* )malloc( width * height * sizeof( double ) );   
	double *dy = ( double* )malloc( width * height * sizeof( double ) );
	double *dxx = ( double* )malloc( width * height * sizeof( double ) );
	double *dyy = ( double* )malloc( width * height * sizeof( double ) );
	double *dxy = ( double* )malloc( width * height  * sizeof( double ) );
	double *Num = ( double* )malloc( width * height  * sizeof( double ) );
	double *Den = ( double* )malloc( width * height  * sizeof( double ) );
	double *I_t = ( double* )malloc( width * height  * sizeof( double ) ); 

	for ( i = 0; i < height; i++)
		for( j = 0; j < width; j++)
			*( tmp + i * width + j ) = *( tmp1 + i * width + j ) = *( I->data + i * width + j ) * 1.0;

	for( m = 0; m < iter; m++){
		//if( flag == 0 ){
			for( i = 0; i < I->height; i++)
				for( j = 0; j < I->width; j++){
					//求dx
					if( j == 0){
						*( dx + i * I->width + j) = ( *( tmp + i * I->width + 1) -  *( tmp + i * I->width ) ) / 2;
					}else if( j == I->width - 1 ){
						*( dx + i * I->width + j) = ( *( tmp + i * I->width + I->width - 1) - *( tmp + i * I->width + I->width - 2) ) / 2;
					}else{
						*( dx + i * I->width + j) = ( *( tmp + i * I->width + ( j + 1 ) ) -  *( tmp + i * I->width + ( j - 1) ) ) / 2;
					}
					//求dy
					if( i == 0 ){
						*( dy + i * I->width + j) = ( *( tmp + 1 * I->width + j ) - *( tmp + j ) ) / 2;
					}else if( i == I->height - 1 ){
						*( dy + i * I->width + j) = ( *( tmp + (I->height - 1) * I->width + j) - *( tmp + (I->height -2) * I->width  + j )) / 2;
					}else{
						*( dy + i * I->width + j) = ( *( tmp + ( i + 1 ) * I->width + j ) - *( tmp + ( i - 1 ) * I->width + j ) ) / 2;
					}
					//求dxx
					if( j == 0 ){
						*( dxx + i * I->width + j) = *( tmp + i * I->width + 1) - *( tmp + i * I->width ) ;
					}else if( j == I->width - 1 ){
						*( dxx + i * I->width + j) = *( tmp + i * I->width + I->width - 2) - *( tmp + i * I->width + I->width - 1 ) ;
					}else{
						*( dxx + i * I->width + j) = *( tmp + i * I->width + ( j + 1 ) ) +  *( tmp + i * I->width + ( j - 1) ) - 2 * ( *( tmp + i * I->width + j));
					}
					//求dyy
					if( i == 0){
						*( dyy + i * I->width + j) = *( tmp + 1 * I->width + j ) - *( tmp + j );
					}else if( i == I->height - 1 ){
						*( dyy + i * I->width + j) = *( tmp + (I->height - 2) * I->width + j) - *( tmp + (I->height -1) * I->width  + j );
					}else{
						*( dyy + i * I->width + j) = *( tmp + ( i + 1 ) * I->width + j ) + *( tmp + ( i - 1 ) * I->width + j ) - 2 *( *( tmp + i * I->width + j ) );
					}
				}//到这里为止，已经把偏导数给求好了
				//这里求dxdy
				for( i = 0; i < I->height; i++)
					for( j = 0; j < I->width; j++){
						if( j == 0 ){
							*( dxy + i * I->width + j ) = ( *( dy + i * I->width + 1 ) - *( dy + i * I->width ) ) / 2;
						}else if( j == I->width - 1 ){
							*( dxy + i * I->width + j ) = ( *( dy + i * I->width + I->width - 1 ) - *( dy + i * I->width + I->width - 2 ) ) / 2;
						}else{
							*( dxy + i * I->width + j ) = ( *( dy + i * I->width + j + 1 ) - *( dy + i * I->width + j -1 ) ) / 2;
						}
				}
				//这里求Num,Den
				for( i = 0; i < I->height; i++)
					for( j = 0; j < I->width; j++){
						*( Num + i * I->width + j ) = ( *( dxx + i * I->width + j ) ) * ( ep + pow( *( dy + i * I->width + j), 2) ) + 
							( *( dyy + i * I->width + j ) ) * ( ep + pow( *( dx + i * I->width + j), 2) ) - 2 * ( *( dx + i * I->width + j)) * ( *( dy + i * I->width + j)) *  ( *( dxy + i * I->width + j));
						*(Den + i * I->width + j ) = pow( (ep + pow( *(dx + i * I->width + j), 2) + pow( *(dy + i * I->width + j), 2) ) , 1.5);
					}

		//这里进行迭代运算
		for( i = 0; i < I->height; i++)
					for( j = 0; j < I->width; j++){
						*( I_t + i * I->width + j ) = (*( Num + i * I->width + j )) / (*( Den + i * I->width + j )) +  0 * (*( tmp1 + i * I->width + j) - *( tmp + i * I->width + j) );
					}
		for( i = 0; i < I->height; i++)
					for( j = 0; j < I->width; j++){
						*( tmp + i * I->width + j ) += dt * (*( I_t + i * I->width + j) );
					}
	}
	for( i = 0; i < I->height; i++)
		for( j = 0; j < I->width; j++){
			*(I->data + i * width + j ) = *( tmp + i * width + j );
		}
	//}
//	flag++;	

		free( dx );
		free( dy );
		free( dxy );
		free( dxx );
		free( dyy );
		free( Num );
		free( Den );
		free( I_t );
		free( tmp );
		free( tmp1 );

}

double calc_lam( Img* I, Img* I0)
{
	int i, j;
	double lam = 0.0, lam1 = 0.0, lam2 = 0.0;
	double *tmp = ( double* )malloc( width * height  * sizeof( double ) );
	double *tmp1 = ( double* )malloc( width * height  * sizeof( double ) );

	double *dx = ( double* )malloc( width * height * sizeof( double ) );   
	double *dy = ( double* )malloc( width * height * sizeof( double ) );
	double *dxx = ( double* )malloc( width * height * sizeof( double ) );
	double *dyy = ( double* )malloc( width * height * sizeof( double ) );
	double *dxy = ( double* )malloc( width * height  * sizeof( double ) );
	double *Num = ( double* )malloc( width * height  * sizeof( double ) );
	double *Den = ( double* )malloc( width * height  * sizeof( double ) );

	for ( i = 0; i < height; i++)
		for( j = 0; j < width; j++)
			*( tmp + i * width + j ) = *( tmp1 + i * width + j ) = *( I->data + i * width + j ) * 1.0;
	
	for( i = 0; i < I->height; i++)
				for( j = 0; j < I->width; j++){
					//求dx
					if( j == 0){
						*( dx + i * I->width + j) = ( *( tmp + i * I->width + 1) -  *( tmp + i * I->width ) ) / 2;
					}else if( j == I->width - 1 ){
						*( dx + i * I->width + j) = ( *( tmp + i * I->width + I->width - 1) - *( tmp + i * I->width + I->width - 2) ) / 2;
					}else{
						*( dx + i * I->width + j) = ( *( tmp + i * I->width + ( j + 1 ) ) -  *( tmp + i * I->width + ( j - 1) ) ) / 2;
					}
					//求dy
					if( i == 0 ){
						*( dy + i * I->width + j) = ( *( tmp + 1 * I->width + j ) - *( tmp + j ) ) / 2;
					}else if( i == I->height - 1 ){
						*( dy + i * I->width + j) = ( *( tmp + (I->height - 1) * I->width + j) - *( tmp + (I->height -2) * I->width  + j )) / 2;
					}else{
						*( dy + i * I->width + j) = ( *( tmp + ( i + 1 ) * I->width + j ) - *( tmp + ( i - 1 ) * I->width + j ) ) / 2;
					}
					//求dxx
					if( j == 0 ){
						*( dxx + i * I->width + j) = *( tmp + i * I->width + 1) - *( tmp + i * I->width ) ;
					}else if( j == I->width - 1 ){
						*( dxx + i * I->width + j) = *( tmp + i * I->width + I->width - 2) - *( tmp + i * I->width + I->width - 1 ) ;
					}else{
						*( dxx + i * I->width + j) = *( tmp + i * I->width + ( j + 1 ) ) +  *( tmp + i * I->width + ( j - 1) ) - 2 * ( *( tmp + i * I->width + j));
					}
					//求dyy
					if( i == 0){
						*( dyy + i * I->width + j) = *( tmp + 1 * I->width + j ) - *( tmp + j );
					}else if( i == I->height - 1 ){
						*( dyy + i * I->width + j) = *( tmp + (I->height - 2) * I->width + j) - *( tmp + (I->height -1) * I->width  + j );
					}else{
						*( dyy + i * I->width + j) = *( tmp + ( i + 1 ) * I->width + j ) + *( tmp + ( i - 1 ) * I->width + j ) - 2 *( *( tmp + i * I->width + j ) );
					}
				}//到这里为止，已经把偏导数给求好了
				//这里求dxdy
				for( i = 0; i < I->height; i++)
					for( j = 0; j < I->width; j++){
						if( j == 0 ){
							*( dxy + i * I->width + j ) = ( *( dy + i * I->width + 1 ) - *( dy + i * I->width ) ) / 2;
						}else if( j == I->width - 1 ){
							*( dxy + i * I->width + j ) = ( *( dy + i * I->width + I->width - 1 ) - *( dy + i * I->width + I->width - 2 ) ) / 2;
						}else{
							*( dxy + i * I->width + j ) = ( *( dy + i * I->width + j + 1 ) - *( dy + i * I->width + j -1 ) ) / 2;
						}
				}
				//这里求Num,Den
				for( i = 0; i < I->height; i++)
					for( j = 0; j < I->width; j++){
						*( Num + i * I->width + j ) = ( *( dxx + i * I->width + j ) ) * ( ep + pow( *( dy + i * I->width + j), 2) ) + 
							( *( dyy + i * I->width + j ) ) * ( ep + pow( *( dx + i * I->width + j), 2) ) - 2 * ( *( dx + i * I->width + j)) * ( *( dy + i * I->width + j)) *  ( *( dxy + i * I->width + j));
						*(Den + i * I->width + j ) = pow( (ep + pow( *(dx + i * I->width + j), 2) + pow( *(dy + i * I->width + j), 2) ) , 1.5);
					}

						//这里计算保真度项系数lam
						for( i = 0; i < I->height; i++)
							for( j = 0; j < I->width; j++)
								lam1 += *( Num + i * I->width + j ) / *( Den + i * I->width + j ) * ( *( tmp + i * I->width + j ) - *( tmp1 + i * I->width + j) );
						lam1 /= width * height;

						for( i = 0; i < I->height; i++)
							for( j = 0; j < I->width; j++)
								lam2 += pow( *( tmp + i * I->width + j ) - *( tmp1 + i * I->width + j ), 2);
						lam2 /= width *height;

						lam = lam1 / lam2;     //第一次计算时，

		free( dx );
		free( dy );
		free( dxy );
		free( dxx );
		free( dyy );
		free( Num );
		free( Den );
		free( tmp );
		free( tmp1 );

		return lam;
}

void tv_demo( Img* I, Img* I0)
{
	double lam;
	double delt;     //计算前后两次迭代结果的误差大小
	Img* old;



	iter = 80;   //迭代次数
	dt = 0.2;   //时间步长
	ep = 1;    //加入的常数值，防止奇异方差

	while( mean( dotSub( I->data, old->data, height, width )) > ){
		
	}
}

	