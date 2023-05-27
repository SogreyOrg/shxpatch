#pragma once
#include "FontLib.h"
class TTF_table;
//CAD 矢量字体 
class fontchrlink
{
public://数据	
	unsigned short  code; //字体编号	
	short  defsz; //数据的长度（def 的长度）	
	char  *def;   //矢量字体信息

	float charWidth;//字宽
	char* symbolName;//符号名
	CharData* m_dataPtr;
public:
	fontchrlink(void);
	~fontchrlink(void);
	//把字体信息转换为矢量字体
	bool ShapeCreateVec(CharData* pOut);
};
