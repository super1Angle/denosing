#include "tv.h"
#include <math.h>


int iter;   //迭代次数
double dt;   //时间步长
int ep;    //加入的常数值，防止奇异方差
double *tmp, *tmp1;
//声明一些偏微分
/***************************************************/
//int flag = 0;
//int i, j;

void tv( double* I, double* I0, double lam)
{
	int i, j, m;

	double *dx = ( double* )malloc( width * height * sizeof( double ) );   
	double *dy = ( double* )malloc( width * height * sizeof( double ) );
	double *dxx = ( double* )malloc( width * height * sizeof( double ) );
	double *dyy = ( double* )malloc( width * height * sizeof( double ) );
	double *dxy = ( double* )malloc( width * height  * sizeof( double ) );
	double *Num = ( double* )malloc( width * height  * sizeof( double ) );
	double *Den = ( double* )malloc( width * height  * sizeof( double ) );
	double *I_t = ( double* )malloc( width * height  * sizeof( double ) ); 

	for( m = 0; m < iter; m++){
			for( i = 0; i < height; i++)
				for( j = 0; j < width; j++){
					//求dx
					if( j == 0){
						*( dx + i * width + j) = ( *( tmp + i * width + 1) -  *( tmp + i * width ) ) / 2;
					}else if( j == width - 1 ){
						*( dx + i * width + j) = ( *( tmp + i * width + width - 1) - *( tmp + i * width + width - 2) ) / 2;
					}else{
						*( dx + i * width + j) = ( *( tmp + i * width + ( j + 1 ) ) -  *( tmp + i * width + ( j - 1) ) ) / 2;
					}
					//求dy
					if( i == 0 ){
						*( dy + i * width + j) = ( *( tmp + 1 * width + j ) - *( tmp + j ) ) / 2;
					}else if( i == height - 1 ){
						*( dy + i * width + j) = ( *( tmp + (height - 1) * width + j) - *( tmp + (height -2) * width  + j )) / 2;
					}else{
						*( dy + i * width + j) = ( *( tmp + ( i + 1 ) * width + j ) - *( tmp + ( i - 1 ) * width + j ) ) / 2;
					}
					//求dxx
					if( j == 0 ){
						*( dxx + i * width + j) = *( tmp + i * width + 1) - *( tmp + i * width ) ;
					}else if( j == width - 1 ){
						*( dxx + i * width + j) = *( tmp + i * width + width - 2) - *( tmp + i * width + width - 1 ) ;
					}else{
						*( dxx + i * width + j) = *( tmp + i * width + ( j + 1 ) ) +  *( tmp + i * width + ( j - 1) ) - 2 * ( *( tmp + i * width + j));
					}
					//求dyy
					if( i == 0){
						*( dyy + i * width + j) = *( tmp + 1 * width + j ) - *( tmp + j );
					}else if( i == height - 1 ){
						*( dyy + i * width + j) = *( tmp + (height - 2) * width + j) - *( tmp + (height -1) * width  + j );
					}else{
						*( dyy + i * width + j) = *( tmp + ( i + 1 ) * width + j ) + *( tmp + ( i - 1 ) * width + j ) - 2 *( *( tmp + i * width + j ) );
					}
				}//到这里为止，已经把偏导数给求好了
				//这里求dxdy
				for( i = 0; i < height; i++)
					for( j = 0; j < width; j++){
						if( j == 0 ){
							*( dxy + i * width + j ) = ( *( dy + i * width + 1 ) - *( dy + i * width ) ) / 2;
						}else if( j == width - 1 ){
							*( dxy + i * width + j ) = ( *( dy + i * width + width - 1 ) - *( dy + i * width + width - 2 ) ) / 2;
						}else{
							*( dxy + i * width + j ) = ( *( dy + i * width + j + 1 ) - *( dy + i * width + j -1 ) ) / 2;
						}
				}
				//这里求Num,Den
				for( i = 0; i < height; i++)
					for( j = 0; j < width; j++){
						*( Num + i * width + j ) = ( *( dxx + i * width + j ) ) * ( ep + pow( *( dy + i * width + j), 2) ) + 
							( *( dyy + i * width + j ) ) * ( ep + pow( *( dx + i * width + j), 2) ) - 2 * ( *( dx + i * width + j)) * ( *( dy + i * width + j)) *  ( *( dxy + i * width + j));
						*(Den + i * width + j ) = pow( (ep + pow( *(dx + i * width + j), 2) + pow( *(dy + i * width + j), 2) ) , 1.5);
					}

		//这里进行迭代运算
		for( i = 0; i < height; i++)
					for( j = 0; j < width; j++){
						*( I_t + i * width + j ) = (*( Num + i * width + j )) / (*( Den + i * width + j )) +  lam * (*( tmp1 + i * width + j) - *( tmp + i * width + j) );
					}
		for( i = 0; i < height; i++)
					for( j = 0; j < width; j++){
						*( tmp + i * width + j ) += dt * (*( I_t + i * width + j) );
					}
	}

		free( dx );
		free( dy );
		free( dxy );
		free( dxx );
		free( dyy );
		free( Num );
		free( Den );
		free( I_t );

}

double calc_lam( double* tmp, double* tmp1)
{
	int i, j;
	double lam = 0.0, lam1 = 0.0, lam2 = 0.0;

	double *dx = ( double* )malloc( width * height * sizeof( double ) );   
	double *dy = ( double* )malloc( width * height * sizeof( double ) );
	double *dxx = ( double* )malloc( width * height * sizeof( double ) );
	double *dyy = ( double* )malloc( width * height * sizeof( double ) );
	double *dxy = ( double* )malloc( width * height  * sizeof( double ) );
	double *Num = ( double* )malloc( width * height  * sizeof( double ) );
	double *Den = ( double* )malloc( width * height  * sizeof( double ) );
	
	for( i = 0; i < height; i++)
				for( j = 0; j < width; j++){
					//求dx
					if( j == 0){
						*( dx + i * width + j) = ( *( tmp + i * width + 1) -  *( tmp + i * width ) ) / 2;
					}else if( j == width - 1 ){
						*( dx + i * width + j) = ( *( tmp + i * width + width - 1) - *( tmp + i * width + width - 2) ) / 2;
					}else{
						*( dx + i * width + j) = ( *( tmp + i * width + ( j + 1 ) ) -  *( tmp + i * width + ( j - 1) ) ) / 2;
					}
					//求dy
					if( i == 0 ){
						*( dy + i * width + j) = ( *( tmp + 1 * width + j ) - *( tmp + j ) ) / 2;
					}else if( i == height - 1 ){
						*( dy + i * width + j) = ( *( tmp + (height - 1) * width + j) - *( tmp + (height -2) * width  + j )) / 2;
					}else{
						*( dy + i * width + j) = ( *( tmp + ( i + 1 ) * width + j ) - *( tmp + ( i - 1 ) * width + j ) ) / 2;
					}
					//求dxx
					if( j == 0 ){
						*( dxx + i * width + j) = *( tmp + i * width + 1) - *( tmp + i * width ) ;
					}else if( j == width - 1 ){
						*( dxx + i * width + j) = *( tmp + i * width + width - 2) - *( tmp + i * width + width - 1 ) ;
					}else{
						*( dxx + i * width + j) = *( tmp + i * width + ( j + 1 ) ) +  *( tmp + i * width + ( j - 1) ) - 2 * ( *( tmp + i * width + j));
					}
					//求dyy
					if( i == 0){
						*( dyy + i * width + j) = *( tmp + 1 * width + j ) - *( tmp + j );
					}else if( i == height - 1 ){
						*( dyy + i * width + j) = *( tmp + (height - 2) * width + j) - *( tmp + (height -1) * width  + j );
					}else{
						*( dyy + i * width + j) = *( tmp + ( i + 1 ) * width + j ) + *( tmp + ( i - 1 ) * width + j ) - 2 *( *( tmp + i * width + j ) );
					}
				}//到这里为止，已经把偏导数给求好了
				//这里求dxdy
				for( i = 0; i < height; i++)
					for( j = 0; j < width; j++){
						if( j == 0 ){
							*( dxy + i * width + j ) = ( *( dy + i * width + 1 ) - *( dy + i * width ) ) / 2;
						}else if( j == width - 1 ){
							*( dxy + i * width + j ) = ( *( dy + i * width + width - 1 ) - *( dy + i * width + width - 2 ) ) / 2;
						}else{
							*( dxy + i * width + j ) = ( *( dy + i * width + j + 1 ) - *( dy + i * width + j -1 ) ) / 2;
						}
				}
				//这里求Num,Den
				for( i = 0; i < height; i++)
					for( j = 0; j < width; j++){
						*( Num + i * width + j ) = ( *( dxx + i * width + j ) ) * ( ep + pow( *( dy + i * width + j), 2) ) + 
							( *( dyy + i * width + j ) ) * ( ep + pow( *( dx + i * width + j), 2) ) - 2 * ( *( dx + i * width + j)) * ( *( dy + i * width + j)) *  ( *( dxy + i * width + j));
						*(Den + i * width + j ) = pow( (ep + pow( *(dx + i * width + j), 2) + pow( *(dy + i * width + j), 2) ) , 1.5);
					}

						//这里计算保真度项系数lam
						for( i = 0; i < height; i++)
							for( j = 0; j < width; j++)
								lam1 += *( Num + i * width + j ) / *( Den + i * width + j ) * ( *( tmp + i * width + j ) - *( tmp1 + i * width + j) );
						lam1 /= width * height;

						for( i = 0; i < height; i++)
							for( j = 0; j < width; j++)
								lam2 += pow( *( tmp + i * width + j ) - *( tmp1 + i * width + j ), 2);
						lam2 /= width *height;

						lam = lam1 / lam2;     //第一次计算时，

		free( dx );
		free( dy );
		free( dxy );
		free( dxx );
		free( dyy );
		free( Num );
		free( Den );

		return lam;
}

void tv_demo( Img* I, Img* I0)
{
	double lam;
	double delt;     //计算前后两次迭代结果的误差大小
	int i, j;
	double* old;
	double threshold = 0.01;

	int  times = 0;      //记录循环执行的次数

	tmp = ( double* )malloc( width * height  * sizeof( double ) );
	tmp1 = ( double* )malloc( width * height  * sizeof( double ) );
	old = ( double* )malloc( width * height  * sizeof( double ) );

	for ( i = 0; i < height; i++)
		for( j = 0; j < width; j++)
			*( tmp + i * width + j ) = *( tmp1 + i * width + j ) = *( I->data + i * width + j ) * 1.0;

	for ( i = 0; i < height; i++)
		for( j = 0; j < width; j++)
			*( old + i * width + j ) = 0.0;

	delt = 1.0;  //这是随便设置的一个值，主要是为了能进入循环
	iter = 80;   //迭代次数
	dt = 0.2;   //时间步长
	ep = 1;    //加入的常数值，防止奇异方差
	lam = 0.0;

	

	while( delt > threshold /*TRUE*/){
		for ( i = 0; i < height; i++)
			for( j = 0; j < width; j++)
				*( old + i * I->width + j )= *( tmp+ i * width + j );
		tv( tmp, tmp1, lam );
		lam = calc_lam( tmp, tmp1 );
		//printf( "%f  \n", lam );
		delt = mean( dotSub( tmp, old, height, width ), height, width );
		//printf( "%f  \n",delt );
		times++;
	}

	printf( "%d", times );

	for( i = 0; i < I->height; i++)
		for( j = 0; j < I->width; j++){
			*(I->data + i * width + j ) = *( tmp + i * width + j );
		}

		free( tmp );
		free( tmp1 );
		free( old );
}

	