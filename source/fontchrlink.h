#pragma once
#include "FontLib.h"
class TTF_table;
//CAD ʸ������ 
class fontchrlink
{
public://����	
	unsigned short  code; //������	
	short  defsz; //���ݵĳ��ȣ�def �ĳ��ȣ�	
	char  *def;   //ʸ��������Ϣ

	float charWidth;//�ֿ�
	char* symbolName;//������
	CharData* m_dataPtr;
public:
	fontchrlink(void);
	~fontchrlink(void);
	//��������Ϣת��Ϊʸ������
	bool ShapeCreateVec(CharData* pOut);
};
