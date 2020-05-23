#pragma once

#define _CRT_SECURE_NO_WARNINGS
#define IN
#define OUT
#define SIZE 8

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

unsigned char getHex(unsigned char ch);
void convertStr2Byte(unsigned char* from, int size, unsigned char* to);

typedef unsigned char byte;
typedef unsigned int word;
typedef unsigned long long ull;

typedef struct POINT {
	word x[2 * SIZE];
	word y[2 * SIZE];
	word z[2 * SIZE];
	int flag;	// 1이면 무한원점, 0이면 값이 있는거
}point;

byte getHex(byte ch);
void convertStr2Byte(byte* from, int size, byte* to);
void readData(word* opA, FILE* fp);
int cmpData(word* A, word* B);

//큰 정수 연산
int addition(IN word* opA, IN word* opB, OUT word* opC);
void additionFp(IN word* opA, IN word* opB, IN word* p, OUT word* opC);
int subtraction(IN word* opA, IN word* opB, OUT word* opC);
void subtractionFp(IN word* opA, IN word* opB, IN word* p, OUT word* opC);
void multiplicationOS(IN word* opA, IN word* opB, OUT word* opC);
void multiplicationPS(IN word* opA, IN word* opB, OUT word* opC);			//PS가 더 빠르다
void mulOSdiv(IN word* opA, IN word* opB, OUT word* opC);
void mulPSdiv(IN word* opA, IN word* opB, OUT word* opC);
void sqr(IN word* opA, OUT word* opC);
void fastReduction(IN word* opC, IN word* p, OUT word* modC);

//역원
void fermatLT(IN word* opZ, IN word* p, OUT word* invZ);
void div2(IN word* opA);
int binaryInversion(IN word* opA, IN word* p, OUT word* opInv);

//타원곡선연산
int ECADD(IN point* P, IN point* Q, IN word* p, OUT point* R);
int ECDBL(IN point* P, IN word* p, OUT point* R);
void LtR(IN word* k, IN point* P, IN word* p, OUT point* kP);
void RtL(IN word* k, IN point* P, IN word* p, OUT point* kP);

//자코비안
int jacobianECDBL(IN point* P, IN word* p, OUT point* R);
int jacobianECADD(IN point* P, IN point* Q, IN word* p, OUT point* R);
void jac2Aff(IN point* P, IN word* p, OUT point* R);
void jacobianLtR(IN word* k, IN point* P, IN word* p, OUT point* R);

//wNAF
void subtractionWord(IN word* opA, IN char opB, OUT word* opC);
void additionWord(word* opA, word opB, word* opC);
void wNAF(IN word* k, IN int w, OUT char* wNAFk);
void makeTable(IN point* P, IN int w, IN word* p, OUT point* table[]);
void wNAFSM(IN point* table[], IN char* wNAFk, IN word* p, OUT point* kP);

//comb method
void precomputationComb(IN point* P, IN int w, IN int num, IN word* p, OUT point* precomputation);
void combPadding(IN word* k, IN int w, OUT char* combK[]);
void combSM(IN word* k, IN point* precomputation, IN word* p, IN char* combK[], IN int w, OUT point* kP);
