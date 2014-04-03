#include "image.h"


Img* read_img( char* filename)
{
	int m_height, m_width, m_component;    //�����洢ͼ����У��У��ͳɷ���
	unsigned char* tmp;
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	JSAMPROW row_pointer[ 1 ];
	Img* image = ( Img* )malloc( sizeof( Img ) );     //imageʹ��֮ǰһ��Ҫ��ʼ��

	FILE* fp = fopen( filename, "rb");
	if( fp == NULL)
		{
		printf( "Can't Open %s", filename );
		exit( 1 );   //�޷����ļ�������������˳�
	}

	cinfo.err = jpeg_std_error( &jerr );
	jpeg_create_decompress( &cinfo );

	jpeg_stdio_src( &cinfo, fp );

	jpeg_read_header( &cinfo, TRUE );

	m_height = cinfo.image_height;
	m_width = cinfo.image_width;
	m_component = cinfo.num_components;
	tmp = ( unsigned char* ) malloc( m_height * m_width * m_component );

	jpeg_start_decompress( &cinfo );


	while( cinfo.output_scanline < cinfo.output_height ){
		row_pointer[ 0 ] = &tmp[ cinfo.output_scanline * m_width * m_component ];
		jpeg_read_scanlines( &cinfo, row_pointer, 1 );
	}

	image->data = tmp;
	image->height = m_height;
	image->width = m_width;
	image->component = m_component;

	height = m_height;
	width = m_width;

	jpeg_finish_decompress( &cinfo );
	jpeg_destroy_decompress( &cinfo );
	fclose( fp );

	return image;
}


Img* rgb2ycrcb( Img* in)
{	
	const float YCbCrYRF = 0.299F;              // RGBתYCbCr��ϵ��(�������ͣ�
    const float YCbCrYGF = 0.587F;
    const float YCbCrYBF = 0.114F;
    const float YCbCrCbRF = -0.168736F;        
    const float YCbCrCbGF = -0.331264F;
    const float YCbCrCbBF = 0.500000F;
    const float YCbCrCrRF = 0.500000F;
    const float YCbCrCrGF = -0.418688F;
    const float YCbCrCrBF = -0.081312F;

	const int Shift = 20;     //�Ŵ���
	const int HalfShiftValue = 1 << (Shift - 1);

    const int YCbCrYRI = (int)(YCbCrYRF * (1 << Shift) + 0.5);         // RGBתYCbCr��ϵ��(�������ͣ�
    const int YCbCrYGI = (int)(YCbCrYGF * (1 << Shift) + 0.5);
    const int YCbCrYBI = (int)(YCbCrYBF * (1 << Shift) + 0.5);
    const int YCbCrCbRI = (int)(YCbCrCbRF * (1 << Shift) + 0.5);
    const int YCbCrCbGI = (int)(YCbCrCbGF * (1 << Shift) + 0.5);
    const int YCbCrCbBI = (int)(YCbCrCbBF * (1 << Shift) + 0.5);
    const int YCbCrCrRI = (int)(YCbCrCrRF * (1 << Shift) + 0.5);
    const int YCbCrCrGI = (int)(YCbCrCrGF * (1 << Shift) + 0.5);
    const int YCbCrCrBI = (int)(YCbCrCrBF * (1 << Shift) + 0.5);


	int i,j;
	int Red, Green, Blue;     //�����洢RGB����
	Img* out = ( Img* )malloc( sizeof( Img ) );
	out->data = ( unsigned char* )malloc( in->height * in->width * 3 );

	out->height =  in->height;
	out->width = in->width;
	out->component = in->component;

	for(i = 0; i < in->height; i++ )
		for(j = 0; j < in->width; j++)
		{
			Red = *( in->data + 3 * ( i * in->width + j ) + 0 );
			Green = *( in->data + 3 * ( i * in->width + j ) + 1 );
			Blue = *( in->data + 3 * ( i * in->width + j ) + 2 );
			*( out->data + 3 * ( i * in->width + j ) + 0 ) = ( unsigned char )(( YCbCrYRI * Red + YCbCrYGI * Green + YCbCrYBI *Blue +HalfShiftValue )>>Shift );
			*( out->data + 3 * ( i * in->width + j ) + 1 ) = ( unsigned char )( 128 +( ( YCbCrCbRI * Red + YCbCrCbGI * Green + YCbCrCbBI *Blue +HalfShiftValue )>>Shift ));
			*( out->data + 3 * ( i * in->width + j ) + 2 ) = ( unsigned char )( 128 +( ( YCbCrCrRI * Red + YCbCrCrGI * Green + YCbCrCrBI *Blue +HalfShiftValue )>>Shift ));
		}

		free( in );    //��Ϊ�Ժ��õ�����YCbCr�ռ䣬���԰�RGB�ռ��ͷŵ�
		return out;
}

Img* get_gray( Img* in)
{
	int i,j;
	Img* out = ( Img* )malloc( sizeof( Img ) );
	out->data = ( unsigned char*)malloc( in->height * in->width );
	for( i = 0; i < in->height; i++)
		for( j = 0; j < in->width; j++)
			*(out->data + i * in->width + j ) = *( in->data + 3 * ( i * in->width + j) );
	out->height = in->height;
	out->width = in->width;
	out->component = 1;      //�Ҷ�ͼ���ͨ������1

	free( in );
	return out;
}


void write_img( char* filename, Img* in)
{
	unsigned char* gray;
	int row_stride;   //ÿһ�е��ֽ���
	FILE* fp;

	JSAMPROW row_pointer[ 1 ];     //һ��λͼ
	struct jpeg_compress_struct jcs;

	struct jpeg_error_mgr jem;
	jcs.err = jpeg_std_error( &jem );
	jpeg_create_compress( &jcs );

	fp = fopen( filename, "wb");
	if( fp == NULL ){
		printf( "Can't Open %s", filename );
		exit( 1 );   //�޷����ļ�������������˳�
	}

	jpeg_stdio_dest( &jcs, fp );

	//����ѹ����������Ҫ������ͼ��Ŀ��ߡ�ɫ��ͨ����ɫ�ʿռ�
	jcs.image_height = in->height;
	jcs.image_width = in->width;
	jcs.input_components = in->component;    //�˴������1��������ǻҶ�ͼ��
	jcs.in_color_space = JCS_GRAYSCALE;

	jpeg_set_defaults( &jcs );
	jpeg_set_quality( &jcs, 80, TRUE );

	jpeg_start_compress( &jcs, TRUE );
	
	row_stride = jcs.image_width;
	gray = in->data;

	while( jcs.next_scanline < jcs.image_height ){
		row_pointer[ 0 ] = &gray[ jcs.next_scanline * row_stride ];
		jpeg_write_scanlines( &jcs, row_pointer, 1 );
	}

	jpeg_finish_compress( &jcs );
	jpeg_destroy_compress( &jcs );
	fclose( fp );
}