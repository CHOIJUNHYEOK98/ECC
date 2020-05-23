#if 1
#include "type.h"

byte getHex(byte ch)
{
	byte hex = 0;

	if (ch >= '0' && ch <= '9')
	{
		hex = ch - '0';
	}
	else if (ch >= 'a' && ch <= 'f')
	{
		hex = ch - 'a' + 10;
	}
	else if (ch >= 'A' && ch <= 'F')
	{
		hex = ch - 'A' + 10;
	}

	return hex;
}
void print256bit(word* A)
{
	int i;
	for (i = SIZE-1; i >= 0; i--)
		printf("%08X", A[i]);
	printf("\n");
}
void convertStr2Byte(byte* from, int size, byte* to)
{
	int cnt_i = 0;
	int cnt_j = 0;
	int ch = 0;

	for (cnt_i = 0; cnt_i < size; cnt_i += 2, cnt_j++)
	{
		ch = from[cnt_i];

		to[cnt_j] = getHex(ch);
		to[cnt_j] <<= 4;

		ch = from[cnt_i + 1];
		to[cnt_j] |= getHex(ch);
	}
}

void readData(word* opA, FILE* fp)
{
	byte buf[100];
	byte buf2[64];
	int i;

	fgets(buf, sizeof(buf), fp);
	convertStr2Byte(buf, 64, buf2);

	for (i = 0; i < SIZE; i++)
		opA[7 - i] = ((buf2[4 * i]) << 24) | ((buf2[4 * i + 1]) << 16) | ((buf2[4 * i + 2]) << 8) | buf2[4 * i + 3];
}

//A와 B를 비교해서 A가 더 크거나 같으면 1을 return하고 아니면 0을 return
int cmpData(word * A, word * B)
{
	int i;
	for (i = SIZE - 1; i >= 0; i--)
	{
		if (A[i] > B[i])
			return 1;
		else if (A[i] < B[i])
			return 0;
	}
	return 1;
}

int addition(IN word* opA, IN word* opB, OUT word* opC)
{
	int i, carry = 0;
	word T = 0;
	for (i = 0; i < SIZE; i++)
	{
		opC[i] = opA[i] + opB[i] + carry;

		if ((carry) && (opC[i] <= opA[i]))
			carry = 1;
		else if (opC[i] < opA[i])
			carry = 1;
		else
			carry = 0;
	}
	return carry;
}

void additionFp(IN word* opA, IN word* opB, IN word* p, OUT word* opC)
{
	if (addition(opA, opB, opC) || cmpData(opC, p))
	{
		subtraction(opC, p, opC);
	}
}

int subtraction(IN word* opA, IN word* opB, OUT word* opC)
{
	int i, borrow = 0;
	for (i = 0; i < SIZE; i++)
	{
		opC[i] = opA[i] - opB[i] - borrow;
		if ((borrow) && (opA[i] <= opB[i]))
			borrow = 1;
		else if (opA[i] < opB[i])
			borrow = 1;
		else
			borrow = 0;
	}
	return borrow;
}

void subtractionFp(IN word* opA, IN word* opB, IN word* p, OUT word* opC)
{
	if (subtraction(opA, opB, opC) || cmpData(opC, p))
		addition(p, opC, opC);
}

//Operand Scanning 방법을 활용한 큰 정수 곱셈
void multiplicationOS(IN word* opA, IN word* opB, OUT word* opC)
{
	ull UV = 0;
	word U = 0;
	memset(opC, 0x00000000, sizeof(word) * (2*SIZE));
	int i, j;	// i : A의 index, j : B의 index
	for (i = 0; i < SIZE; i++)
	{
		U = 0x00000000;
		for (j = 0; j < SIZE; j++)
		{
			UV = (ull)opC[i + j] + ((ull)opA[i] * (ull)opB[j]) + (ull)U;

			opC[i + j] = (UV & 0xffffffff);
			U = (word)((UV >> 32));
		}
		opC[i + SIZE] = U;
	}
}
#if 1
void multiplicationPS(IN word* opA, IN word* opB, OUT word* opC)
{
	word R0 = 0, R1 = 0, R2 = 0, U = 0, V = 0, temp = 0;
	ull UV = 0;
	int i, j, k, carry = 0;
	for (k = 0; k < (2*SIZE-1); k++)
	{
		i = 0; j = k;
		while (j > (SIZE-1))
		{
			i++; j--;
		}
		while ((j >= 0) && (j <= (SIZE-1)) && (i >= 0) && (i <= (SIZE-1)))
		{
			carry = 0;

			UV = (ull)opA[i] * (ull)opB[j];
			V = (word)((UV & 0xffffffff));
			U = (word)((UV >> 32));

			R0 += V;
			if (R0 < V)
				carry = 1;
			else
				carry = 0;

			temp = R1 + U;

			if (temp < R1)
			{
				R1 = temp + carry;
				carry = 1;
			}
			else if ((temp == 0xffffffff) && (carry))
			{
				R1 = temp + carry;
				carry = 1;
			}
			else
			{
				R1 = temp + carry;
				carry = 0;
			}

			R2 = R2 + carry;

			i++;
			j--;
		}
		opC[k] = R0; R0 = R1; R1 = R2; R2 = 0;
	}
	opC[2 * SIZE - 1] = R0;
}
void mulOSdiv(IN word* opA, IN word* opB, OUT word* opC) //32bit word를 16bit word로 나누어 곱셈
{
	word UV = 0;	//32bit
	word U = 0;		//16bit

	word A = 0;		//16

	memset(opC, 0x00000000, sizeof(word) * (2 * SIZE));
	int i, j, carry = 0;
	for (i = 0; i < SIZE; i++)	//A[i]의 index
	{
		U = 0x0000;
		A = opA[i] & 0xffff;
		for (j = 0; j < SIZE; j++)	//B[i]의 index
		{
			UV = (opC[i + j] &0xffff) + A * (opB[j] & 0xffff) + U;
			opC[i + j] = (opC[i + j] & 0xffff0000) + (UV & 0xffff);
			U = UV >> 16;

			UV = (opC[i + j] >> 16) + A * (opB[j] >> 16) + U;
			opC[i + j] = (opC[i + j] & 0xffff) + ((UV & 0xffff) << 16);
			U = UV >> 16;
		}
		opC[i + SIZE] = (opC[i + SIZE] & 0xffff0000) + U;

		U = 0x0000;
		A = opA[i] >> 16;
		for (j = 0; j < SIZE; j++)	//B[i]의 index
		{
			UV = (opC[i + j] >> 16) + A * (opB[j] & 0xffff) + U;
			opC[i + j] = (opC[i + j] & 0xffff) + ((UV & 0xffff) << 16);
			U = UV >> 16;
			
			UV = (opC[i + j + 1] & 0xffff) + A * (opB[j] >> 16) + U;
			opC[i + j + 1] = (opC[i + j + 1] & 0xffff0000) + (UV & 0xffff);
			U = UV >> 16;
		}
		opC[i + SIZE] = (opC[i + SIZE] & 0xffff) + (U << 16);
	}
}
void mulPSdiv(IN word* opA, IN word* opB, OUT word* opC)
{
	word R0 = 0, R1 = 0, R2 = 0, R3 = 0, R4 = 0;
	word U = 0, V = 0, temp = 0, A = 0, B = 0;	//16bit
	word UV = 0; //32bit
	int i, j, k, carry = 0;
	for (k = 0; k < (2 * SIZE - 1); k++)
	{
		i = 0; j = k;
		while (j > (SIZE - 1))
		{
			i++; j--;
		}
		while ((j >= 0) && (j <= (SIZE - 1)) && (i >= 0) && (i <= (SIZE - 1)))
		{
			carry = 0;
			// A[i], B[j]의 하위 16비트의 곱
			A = opA[i] & 0xffff;
			B = opB[j] & 0xffff;

			UV = A * B;
			V = UV & 0xffff;
			U = UV >> 16;

			R0 = (R0 + V) & 0xffff;
			if (R0 < V)
				carry = 1;
			else
				carry = 0;

			temp = (R1 + U) & 0xffff;

			if (temp < R1)
			{
				R1 = temp + carry;
				carry = 1;
			}
			else if ((temp == 0xffff) && (carry))
			{
				R1 = (temp + carry) & 0xffff;
				carry = 1;
			}
			else
			{
				R1 = temp + carry;
				carry = 0;
			}

			R2 = (R2 + carry) & 0xffff;

			//A[i]의 하위 16비트와 B[j]의 상위 16비트의 곱
			B = opB[j] >> 16;

			UV = A * B;
			V = UV & 0xffff;
			U = UV >> 16;

			R1 = (R1 + V) & 0xffff;
			if (R1 < V)
				carry = 1;
			else
				carry = 0;

			temp = (R2 + U) & 0xffff;

			if (temp < R2)
			{
				R2 = temp + carry;
				carry = 1;
			}
			else if ((temp == 0xffff) && (carry))
			{
				R2 = (temp + carry) & 0xffff;
				carry = 1;
			}
			else
			{
				R2 = temp + carry;
				carry = 0;
			}

			R3 = (R3 + carry) & 0xffff;

			//A[i]의 상위 16비트와 B[j]의 하위 16비트의 곱
			A = opA[i] >> 16;
			B = opB[j] & 0xffff;

			UV = A * B;
			V = UV & 0xffff;
			U = UV >> 16;

			R1 = (R1 + V) & 0xffff;
			if (R1 < V)
				carry = 1;
			else
				carry = 0;

			temp = (R2 + U) & 0xffff;

			if (temp < R2)
			{
				R2 = temp + carry;
				carry = 1;
			}
			else if ((temp == 0xffff) && (carry))
			{
				R2 = (temp + carry) & 0xffff;
				carry = 1;
			}
			else
			{
				R2 = temp + carry;
				carry = 0;
			}

			R3 = (R3 + carry) & 0xffff;

			//A[i]의 상위 16비트와 B[j]의 상위 16비트의 곱
			B = opB[j] >> 16;

			UV = A * B;
			V = UV & 0xffff;
			U = UV >> 16;

			R2 = (R2 + V) & 0xffff;
			if (R2 < V)
				carry = 1;
			else
				carry = 0;

			temp = (R3 + U) & 0xffff;

			if (temp < R3)
			{
				R3 = temp + carry;
				carry = 1;
			}
			else if ((temp == 0xffff) && (carry))
			{
				R3 = (temp + carry) & 0xffff;
				carry = 1;
			}
			else
			{
				R3 = temp + carry;
				carry = 0;
			}

			R4 = (R4 + carry) & 0xffff;

			i++;
			j--;
		}
		opC[k] = ((R1 & 0xffff) << 16) + (R0 & 0xffff);
		R0 = R2; R1 = R3; R2 = R4; R3 = 0; R4 = 0;
	}
	opC[2 * SIZE - 1] = ((R1 & 0xffff) << 16) + (R0 & 0xffff);
}
void sqr(IN word* opA, OUT word* opC)	//구현 진행 중
{
	word R0 = 0, R1 = 0, R2 = 0, U = 0, V = 0, temp = 0;
	ull UV = 0;
	int i, j, k, carry = 0;
	for (k = 0; k < (2 * SIZE - 1); k++)
	{
		i = 0; j = k;
		while (j > (SIZE - 1))
		{
			i++; j--;
		}
		while ((j >= 0) && (j <= (SIZE - 1)) && (i >= 0) && (i <= (SIZE - 1)))
		{
			UV = (ull)opA[i] * (ull)opA[j];
			if (i < j)
			{
				if (UV >> 63)
					R2++;
				UV <<= 1;
			}
			else if(i > j)
			{
				break;
			}
			V = (word)((UV & 0xffffffff));
			U = (word)((UV >> 32));

			R0 += V;
			if (R0 < V)
				carry = 1;
			else
				carry = 0;

			temp = R1 + U;

			if (temp < R1)
			{
				R1 = temp + carry;
				carry = 1;
			}
			else if ((temp == 0xffffffff) && (carry))
			{
				R1 = temp + carry;
				carry = 1;
			}
			else
			{
				R1 = temp + carry;
				carry = 0;
			}

			R2 = R2 + carry;

			i++;
			j--;
		}
		opC[k] = R0; R0 = R1; R1 = R2; R2 = 0;
	}
	opC[2 * SIZE - 1] = R0;
}
void fastReduction(IN word* opC, IN word* p, OUT word* modC)
{
	word s[9][8] = { {0,}, };
	word temp[2][8] = { 0, };
	word out[8] = { 0, };
	s[0][0] = opC[0]; s[0][1] = opC[1]; s[0][2] = opC[2]; s[0][3] = opC[3];
	s[0][4] = opC[4]; s[0][5] = opC[5]; s[0][6] = opC[6]; s[0][7] = opC[7];

	s[1][3] = opC[11]; s[1][4] = opC[12]; s[1][5] = opC[13];
	s[1][6] = opC[14]; s[1][7] = opC[15];

	s[2][3] = opC[12]; s[2][4] = opC[13];
	s[2][5] = opC[14]; s[2][6] = opC[15];

	s[3][0] = opC[8]; s[3][1] = opC[9]; s[3][2] = opC[10];
	s[3][6] = opC[14]; s[3][7] = opC[15];

	s[4][0] = opC[9]; s[4][1] = opC[10]; s[4][2] = opC[11]; s[4][3] = opC[13];
	s[4][4] = opC[14]; s[4][5] = opC[15]; s[4][6] = opC[13]; s[4][7] = opC[8];

	s[5][0] = opC[11]; s[5][1] = opC[12]; s[5][2] = opC[13];
	s[5][6] = opC[8]; s[5][7] = opC[10];

	s[6][0] = opC[12]; s[6][1] = opC[13]; s[6][2] = opC[14];
	s[6][3] = opC[15]; s[6][6] = opC[9]; s[6][7] = opC[11];

	s[7][0] = opC[13]; s[7][1] = opC[14]; s[7][2] = opC[15]; s[7][3] = opC[8];
	s[7][4] = opC[9]; s[7][5] = opC[10]; s[7][7] = opC[12];

	s[8][0] = opC[14]; s[8][1] = opC[15]; s[8][3] = opC[9];
	s[8][4] = opC[10]; s[8][5] = opC[11]; s[8][7] = opC[13];

	additionFp(s[1], s[1], p, temp[0]);
	additionFp(s[2], s[2], p, temp[1]);
	additionFp(s[0], temp[0], p, out);
	additionFp(out, temp[1], p, s[0]);
	additionFp(s[0], s[3], p, s[1]);
	additionFp(s[1], s[4], p, s[2]);
	subtractionFp(s[2], s[5], p, s[3]);
	subtractionFp(s[3], s[6], p, s[4]);
	subtractionFp(s[4], s[7], p, s[5]);
	subtractionFp(s[5], s[8], p, modC);
}
void fermatLT(IN word* opZ, IN word* p, OUT word* invZ)
{
	word z3[2 * SIZE] = { 0, }, t0[2 * SIZE] = { 0, }, t1[2 * SIZE] = { 0, },
		t2[2 * SIZE] = { 0, }, t3[2 * SIZE] = { 0, }, temp[2 * SIZE] = { 0, }, temp2[2 * SIZE] = { 0, };
	int i;
	multiplicationPS(opZ, opZ, temp);
	fastReduction(temp, p, temp2);
	multiplicationPS(temp2, opZ, temp);
	fastReduction(temp, p, z3);

	multiplicationPS(z3, z3, temp);
	fastReduction(temp, p, temp2);
	multiplicationPS(temp2, temp2, temp);
	fastReduction(temp, p, temp2);
	multiplicationPS(temp2, z3, temp);
	fastReduction(temp, p, temp2);
	//temp2 == z15

	multiplicationPS(temp2, temp2, temp);
	fastReduction(temp, p, temp2);
	multiplicationPS(temp2, temp2, temp);
	fastReduction(temp, p, temp2);
	multiplicationPS(temp2, z3, temp);
	fastReduction(temp, p, t0);

	multiplicationPS(t0, t0, temp);
	fastReduction(temp, p, temp2);
	for (i = 0; i < 5; i++) {
		multiplicationPS(temp2, temp2, temp);
		fastReduction(temp, p, temp2);
	}
	multiplicationPS(temp2, t0, temp);
	fastReduction(temp, p, temp2);
	memcpy(t1, temp2, sizeof(word) * SIZE);
	//temp2 == t1
	
	for (i = 0; i < 12; i++) {
		multiplicationPS(temp2, temp2, temp);
		fastReduction(temp, p, temp2);
	}
	multiplicationPS(temp2, t1, temp);
	fastReduction(temp, p, temp2);
	for (i = 0; i < 6; i++) {
		multiplicationPS(temp2, temp2, temp);
		fastReduction(temp, p, temp2);
	}
	multiplicationPS(temp2, t0, temp);
	fastReduction(temp, p, t2);
	//temp2 == t2

	multiplicationPS(t2, t2, temp);
	fastReduction(temp, p, temp2);
	multiplicationPS(temp2, temp2, temp);
	fastReduction(temp, p, temp2);
	multiplicationPS(temp2, z3, temp);
	fastReduction(temp, p, temp2);
	memcpy(t3, temp2, sizeof(word) * SIZE);

	for (i = 0; i < 32; i++) {
		multiplicationPS(temp2, temp2, temp);
		fastReduction(temp, p, temp2);
	}
	multiplicationPS(temp2, opZ, temp);
	fastReduction(temp, p, temp2);
	for (i = 0; i < 96; i++) {
		multiplicationPS(temp2, temp2, temp);
		fastReduction(temp, p, temp2);
	}
	//temp2 == t4

	for (i = 0; i < 32; i++) {
		multiplicationPS(temp2, temp2, temp);
		fastReduction(temp, p, temp2);
	}
	multiplicationPS(temp2, t3, temp);
	fastReduction(temp, p, temp2);
	for (i = 0; i < 32; i++) {
		multiplicationPS(temp2, temp2, temp);
		fastReduction(temp, p, temp2);
	}
	multiplicationPS(temp2, t3, temp);
	fastReduction(temp, p, temp2);
	//temp2 == t5

	for (i = 0; i < 30; i++) {
		multiplicationPS(temp2, temp2, temp);
		fastReduction(temp, p, temp2);
	}
	multiplicationPS(temp2, t2, temp);
	fastReduction(temp, p, temp2);
	multiplicationPS(temp2, temp2, temp);
	fastReduction(temp, p, temp2);
	multiplicationPS(temp2, temp2, temp);
	fastReduction(temp, p, temp2);
	multiplicationPS(temp2, opZ, temp);
	fastReduction(temp, p, invZ);
}
void div2(IN word* opA)
{
	int i;
	for (i = 0; i < SIZE - 1; i++)
		opA[i] = (((opA[i + 1] & 0x1) << 31) | (opA[i] >> 1));
	opA[i] >>= 1;
}
int binaryInversion(IN word* opA, IN word* p, OUT word* opInv)
{
	word u[SIZE] = { 0, }, v[SIZE] = { 0, }, x1[SIZE] = { 0, }, x2[SIZE] = { 0, },
		cmp[SIZE] = { 0, }, temp[SIZE] = { 0, };
	if (!(memcmp(opA, u, sizeof(word) * SIZE)))
	{
		memset(opInv, 0, sizeof(word) * SIZE);
		return 0;
	}
	int carry = 0;
	memcpy(u, opA, sizeof(word) * SIZE);
	memcpy(v, p, sizeof(word) * SIZE);
	x1[0] = 1; cmp[0] = 1;
	
	while (memcmp(u, cmp, sizeof(word) * SIZE) && memcmp(v, cmp, sizeof(word) * SIZE))
	{
		while (!(u[0] & 0x1))
		{
			div2(u);
			if (!(x1[0] & 0x1))
				div2(x1);
			else
			{
				carry = addition(x1, p, temp);
				div2(temp);
				memcpy(x1, temp, sizeof(word) * SIZE);
				x1[SIZE - 1] += (carry << 31);
			}
		}
		while (!(v[0] & 0x1))
		{
			div2(v);
			if (!(x2[0] & 0x1))
				div2(x2);
			else
			{
				carry = addition(x2, p, temp);
				div2(temp);
				memcpy(x2, temp, sizeof(word) * SIZE);
				x2[SIZE - 1] += (carry << 31);
			}
		}
		if (cmpData(u, v))
		{
			subtractionFp(u, v, p, temp);
			memcpy(u, temp, sizeof(word) * SIZE);
			subtractionFp(x1, x2, p, temp);
			memcpy(x1, temp, sizeof(word) * SIZE);
		}
		else
		{
			subtractionFp(v, u, p, temp);
			memcpy(v, temp, sizeof(word) * SIZE);
			subtractionFp(x2, x1, p, temp);
			memcpy(x2, temp, sizeof(word) * SIZE);
		}
	}
	if (!memcmp(u, cmp, sizeof(word) * SIZE))
		memcpy(opInv, x1, sizeof(word) * SIZE);
	else
		memcpy(opInv, x2, sizeof(word) * SIZE);
	return 0;
}
//return == 0이면 값이 바뀐 것이고, return == 1이면 P, Q 둘 다 무한원점
int ECADD(IN point* P, IN point* Q, IN word* p, OUT point* R)
{
	if ((P->flag == 1) && (Q->flag == 0))
	{
		memcpy(R, Q, sizeof(point));
		return 0;
	}
	else if ((P->flag == 0) && (Q->flag == 1))
	{
		memcpy(R, P, sizeof(point));
		return 0;
	}
	else if ((P->flag == 1) && (Q->flag == 1))
		return 1;

	word temp1[2 * SIZE] = { 0, }, temp2[2 * SIZE] = { 0, }, tempI[2 * SIZE] = { 0, }, tempG[2 * SIZE] = { 0, }, tempM[2 * SIZE] = { 0, };
	word Rx[SIZE] = { 0, }, Ry[SIZE] = { 0, };
	subtractionFp(Q->y, P->y, p, temp1);	//temp1 = y2-y1
	subtractionFp(Q->x, P->x, p, temp2);	//temp2 = x2-x1
	binaryInversion(temp2, p, tempI);		//tempI = 1/(x2-x1)
	multiplicationPS(temp1, tempI, temp2);	//temp2 = (y2-y1)/(x2-x1)
	fastReduction(temp2, p, tempG);			//tempG = (y2-y1)/(x2-x1)의 감산 후 결과
	multiplicationPS(tempG, tempG, tempM);	//tempM = ((y2-y1)/(x2-x1))^2
	fastReduction(tempM, p, temp1);			//temp1 = ((y2-y1)/(x2-x1))^2 의 감산 후 결과
	subtractionFp(temp1, P->x, p, temp2);
	subtractionFp(temp2, Q->x, p, Rx);

	subtractionFp(P->x, Rx, p, temp1);	//temp1 = x1-x3
	multiplicationPS(tempG, temp1, temp2);	//temp2 = tempG*(x1-x3)
	fastReduction(temp2, p, temp1);			//temp1 = temp2의 감산
	subtractionFp(temp1, P->y, p, R->y);

	memcpy(R->x, Rx, sizeof(word) * SIZE);

	return 0;
}
// a = -3
// P가 무한원점이면 return 1
int ECDBL(IN point* P, IN word* p, OUT point* R)
{
	if (P->flag == 1)
		return 1;

	word temp1[2 * SIZE] = { 0, }, temp2[2 * SIZE] = { 0, }, temp3[2 * SIZE] = { 0, }, tempG[2*SIZE] = { 0, }, a[2*SIZE] = { 0, };
	word Rx[SIZE] = { 0, }, Ry[SIZE] = { 0, };
	multiplicationPS(P->x, P->x, temp1);	//temp1 = x1^2
	fastReduction(temp1, p, temp2);			//temp2 = x1^2의 감산
	additionFp(temp2, temp2, p, temp3);		//temp3 = 2*x1^2
	additionFp(temp3, temp2, p, temp1);		//temp1 = 3*x1^2
	memset(temp2, 0, sizeof(word) * SIZE);
	temp2[0] = 3;
	subtractionFp(p, temp2, p, a);			//a = -3 = p-3
	additionFp(temp1, a, p, temp2);			//temp2 = 3*x1^2+a
	additionFp(P->y, P->y, p, temp1);		//temp1 = 2*y1
	binaryInversion(temp1, p, temp3);		//temp3 = 1/(2*y1)
	multiplicationPS(temp2, temp3, temp1);	//temp1 = ((3*x1^2+a)/(2*y1))
	fastReduction(temp1, p, tempG);
	multiplicationPS(tempG, tempG, temp1);
	fastReduction(temp1, p, temp2);			//temp2 = ((3*x1^2+a)/(2*y1))^2의 감산
	additionFp(P->x, P->x, p,temp1);
	subtractionFp(temp2, temp1, p, Rx);

	subtractionFp(P->x, Rx, p, temp1);
	multiplicationPS(tempG, temp1, temp2);
	fastReduction(temp2, p, temp1);
	subtractionFp(temp1, P->y, p, R->y);

	memcpy(R->x, Rx, sizeof(word) * SIZE);

	return 0;
}
//Q = k*P
void LtR(IN word* k, IN point* P, IN word* p, OUT point* kP)
{
	point Q;
	Q.flag = 1;
	int i, j;
	for (i = SIZE - 1; i >= 0; i--)
	{
		for (j = 31; j >= 0; j--)
		{
			Q.flag = ECDBL(&Q, p, &Q);
			if ((k[i] >> j) & 0x1)
			{
				Q.flag = ECADD(&Q, P, p, &Q);
			}
		}
	}
	memcpy(kP, &Q, sizeof(point));
}
void RtL(IN word* k, IN point* P, IN word* p, OUT point* kP)
{
	point Q, temp;
	memcpy(&temp, P, sizeof(point));
	Q.flag = 1;
	int i, j;
	for (i = 0; i < SIZE; i++)
	{
		for (j = 0; j < 32; j++)
		{
			if ((k[i] >> j) & 0x1)
				Q.flag = ECADD(&Q, &temp, p, &Q);
			temp.flag = ECDBL(&temp, p, &temp);
		}
	}
	memcpy(kP, &Q, sizeof(point));
}
int jacobianECDBL(IN point* P, IN word* p, OUT point* R)
{
	if (P->flag)							//1
		return 1;

	word T1[2 * SIZE] = { 0, }, T2[2 * SIZE] = { 0, }, T3[2 * SIZE] = { 0, };
	word temp1[2 * SIZE] = { 0, }, temp2[2 * SIZE] = { 0, }, temp3[2 * SIZE] = { 0, };
	word Rx[SIZE] = { 0, }, Ry[SIZE] = { 0, }, Rz[SIZE] = { 0, };
	int carry = 0;

	multiplicationPS(P->z, P->z, temp1);	//2
	fastReduction(temp1, p, T1);	//T1 = Z1^2

	subtractionFp(P->x, T1, p, T2);			//3

	additionFp(P->x, T1, p, temp1);			//4
	memcpy(T1, temp1, sizeof(word) * SIZE);

	multiplicationPS(T1, T2, temp1);		//5
	fastReduction(temp1, p, T2);

	additionFp(T2, T2, p, temp1);			//6
	additionFp(temp1, T2, p, temp2);
	memcpy(T2, temp2, sizeof(word) * SIZE);

	additionFp(P->y, P->y, p, temp1);		//7, temp1 = Y3

	multiplicationPS(temp1, P->z, temp2);	//8
	fastReduction(temp2, p, R->z);

	multiplicationPS(temp1, temp1, temp2);	//9
	fastReduction(temp2, p, temp1);

	multiplicationPS(temp1, P->x, temp2);	//10
	fastReduction(temp2, p, T3);

	multiplicationPS(temp1, temp1, temp2);	//11
	fastReduction(temp2, p, temp1);

	if (!(temp1[0] & 0x1))						//12
		div2(temp1);
	else
	{
		carry = addition(temp1, p, temp2);
		div2(temp2);
		memcpy(temp1, temp2, sizeof(temp1));
		temp1[SIZE - 1] += (carry << 31);
	}

	multiplicationPS(T2, T2, temp2);		//13
	fastReduction(temp2, p, temp3);	//temp3 = X3

	additionFp(T3, T3, p, temp2);			//14

	subtractionFp(temp3, temp2, p, R->x);	//15

	subtractionFp(T3, R->x, p, temp3);		//16

	multiplicationPS(temp3, T2, temp2);		//17
	fastReduction(temp2, p, T1);

	subtractionFp(T1, temp1, p, R->y);		//18

	return 0;
}
int jacobianECADD(IN point* P, IN point* Q, IN word* p, OUT point* R)
{
	if (Q->flag)
	{
		memcpy(R, P, sizeof(point));
		return 0;
	}
	else if (P->flag)
	{
		memcpy(R, Q, sizeof(point));
		memset(R->z, 0, sizeof(R->z));
		R->z[0] = 1;
		return 0;
	}
	word temp1[2 * SIZE] = { 0, }, temp2[2 * SIZE] = { 0, }, temp3[2 * SIZE] = { 0, }, temp4[2 * SIZE] = { 0, };
	word T1[2 * SIZE] = { 0, }, T2[2 * SIZE] = { 0, }, T3[2 * SIZE] = { 0, }, T4[2 * SIZE] = { 0, };
	word Rx[SIZE] = { 0, }, Ry[SIZE] = { 0, }, Rz[SIZE] = { 0, };

	multiplicationPS(P->z, P->z, temp1);	//3
	fastReduction(temp1, p, T1);			
	multiplicationPS(T1, P->z, temp2);		//4
	fastReduction(temp2, p, T2);			
	multiplicationPS(T1, Q->x, temp1);		//5
	fastReduction(temp1, p, T1);			
	multiplicationPS(T2, Q->y, temp1);		//6
	fastReduction(temp1, p, T2);
	subtractionFp(T1, P->x, p, temp1);		//7
	memcpy(T1, temp1, sizeof(word) * (2 * SIZE));
	subtractionFp(T2, P->y, p, temp1);		//8
	memcpy(T2, temp1, sizeof(word) * (2 * SIZE));
	if (!(memcmp(T1, temp3, sizeof(word) * SIZE)))	//9
	{
		if (!(memcmp(T2, temp3, sizeof(word) * SIZE)))
		{
			memset(Q->z, 0, sizeof(Q->z));
			Q->z[0] = 1;
			jacobianECDBL(Q, p, R);
			return 0;
		}
		else
			return 1;
	}
	multiplicationPS(P->z, T1, temp1);		//10
	fastReduction(temp1, p, R->z);
	multiplicationPS(T1, T1, temp1);		//11
	fastReduction(temp1, p, T3);
	multiplicationPS(T3, T1, temp1);		//12
	fastReduction(temp1, p, T4);
	multiplicationPS(T3, P->x, temp1);		//13
	fastReduction(temp1, p, T3);
	additionFp(T3, T3, p, T1);				//14
	multiplicationPS(T2, T2, temp1);		//15
	fastReduction(temp1, p, Rx);
	subtractionFp(Rx, T1, p, temp1);		//16
	subtractionFp(temp1, T4, p, Rx);		//17
	subtractionFp(T3, Rx, p, temp1);		//18
	multiplicationPS(temp1, T2, temp2);		//19
	fastReduction(temp2, p, T3);
	multiplicationPS(T4, P->y, temp1);		//20
	fastReduction(temp1, p, T4);
	subtractionFp(T3, T4, p, R->y);			//21

	memcpy(R->x, Rx, sizeof(Rx));

	return 0;
}
void jac2Aff(IN point* P, IN word* p, OUT point* R)
{
	word temp1[2 * SIZE] = { 0, }, temp2[2 * SIZE] = { 0, }, tempI[2 * SIZE] = { 0, };
	multiplicationPS(P->z, P->z, temp1);
	fastReduction(temp1, p, temp2);
	multiplicationPS(temp2, P->z, temp1);
	fastReduction(temp1, p, temp2);
	binaryInversion(temp2, p, tempI);
	multiplicationPS(P->y, tempI, temp1);
	fastReduction(temp1, p, R->y);
	multiplicationPS(P->z, tempI, temp1);
	fastReduction(temp1, p, temp2);
	multiplicationPS(P->x, temp2, temp1);
	fastReduction(temp1, p, R->x);
}
void jacobianLtR(IN word* k, IN point* P, IN word* p, OUT point* R)
{
	point Q;
	Q.flag = 1;
	int i, j;
	for (i = SIZE - 1; i >= 0; i--)
	{
		for (j = 31; j >= 0; j--)
		{
			Q.flag = jacobianECDBL(&Q, p, &Q);
			if ((k[i] >> j) & 0x1)
			{
				Q.flag = jacobianECADD(&Q, P, p, &Q);
			}
		}
	}
	jac2Aff(&Q, p, R);
}
//큰 정수 opA, int opB, 결과를 저장 할 opC, opA의 크기 size
void subtractionWord(IN word* opA, IN char opB, OUT word* opC)
{
	int i = 0, borrow = 0;
	word temp[SIZE] = { 0, };
	temp[0] = opA[0] - opB;
	for (i = 1; i < SIZE; i++)
		temp[i] = opA[i];
	memcpy(opC, temp, sizeof(word) * SIZE);
}
void additionWord(word* opA, word opB, word* opC)
{
	int i = 0, carry = 0;
	word temp[SIZE] = { 0, };
	temp[i] = opA[i] + opB;
	if (temp[i] < opA[i])
		carry = 1;
	for (i = 1; i < SIZE; i++)
	{
		temp[i] = opA[i] + carry;
		if (opA[i] != 0xffffffff)
			carry = 0;
	}
	memcpy(opC, temp, sizeof(word) * SIZE);
}
void wNAF(IN word* k, IN int w, OUT char* wNAFk)
{
	int i = 0;
	word modAnd = (word)pow(2, w) - 1;		//mod 2^w를 하기 위한 변수
	word modValue = ((modAnd + 1) >> 1);	//wNAFk값이 +인지 -인지를 결정해주기 위한 2^(w-1)를 저장하는 변수
	word modSub = modValue << 1;			//2^w가 저장되어 있는 변수
	while (k[0] != 0)
	{
		if (k[0] & 0x1)
		{
			wNAFk[i] = (k[0] & modAnd);
			if (wNAFk[i] > modValue)
			{
				wNAFk[i] -= modSub;
			}
			subtractionWord(k, wNAFk[i], k);
		}
		else
			wNAFk[i] = 0;
		div2(k);
		i++;
	}
	for (; i < 257; i++)
		wNAFk[i] = 0;
}
//2^w로 wNAFk를 했을때 나올 수 있는 경우에 대해서 사전계산 table을 만들기 위한 함수
//table[i] = (2*i + 1)P., 총 2^(w-2)개를 저장해야 한다.
void makeTable(IN point* P, IN int w, IN word* p, OUT point* table[])
{
	int i;
	int num = (int)pow(2, (double)w - 2);	//총 2*(w-2)개를 저장해야한다.
	point P2;
	point temp;
	P2.flag = 1;
	P2.flag = ECDBL(P, p, &P2);				//2P값을 P2에 저장한다.
	memcpy(table[0], P, sizeof(point));		//table[0]에 P를 저장한다.

	for (i = 1; i < num; i++)
	{
		table[i]->flag = ECADD(table[i - 1], &P2, p, table[i]);
	}
}
void wNAFSM(IN point* table[], IN char* wNAFk, IN word* p, OUT point* kP)
{
	point Q;
	int i;
	Q.flag = 1;
	point temp;

	for (i = 256; i >= 0; i--)
	{
		Q.flag = jacobianECDBL(&Q, p, &Q);
		//Q.flag = ECDBL(&Q, p, &Q);
		if (wNAFk[i] != 0)
		{
			if (wNAFk[i] > 0)
			{
				Q.flag = jacobianECADD(&Q, table[((wNAFk[i] - 1) / 2)], p, &Q);
				//Q.flag = ECADD(&Q, table[((wNAFk[i] - 1) / 2)], p, &Q);
			}
			else //wNAFk[i] < 0인 경우
			{
				subtractionFp(p, table[(((-wNAFk[i]) - 1) / 2)]->y, p, &temp.y);
				memcpy(&temp.x, table[(((-wNAFk[i]) - 1) / 2)]->x, sizeof(word) * SIZE);
				//Q.flag = ECADD(&Q, &temp, p, &Q);
				Q.flag = jacobianECADD(&Q, &temp, p, &Q);
			}
		}
	}
	jac2Aff(&Q, p, kP);
	//memcpy(kP, &Q, sizeof(point));
}
//2^0, 2^d, 2^(2d), ... 의 가능한 모든 덧셈 조합 2^w가지를 precomputation에 저장 ( d = 256/w )
//num = 2^w
void precomputationComb(IN point* P, IN int w, IN int num, IN word* p, OUT point* precomputation)
{
	int i, j, d, cnt;
	word** temp;			//2^0, 2^d, 2^(2d), ... 를 temp에 저장하고
	point* temp2;			//(2^0)*P, (2^d)*P, (2^(2d))*P, ...를 temp2에 저장
	point Q;
	Q.flag = 1;
	temp = (word**)malloc(sizeof(word*) * w);
	for (i = 0; i < w; i++)
	{
		temp[i] = (word*)malloc(sizeof(word) * SIZE);
		memset(temp[i], 0, sizeof(word) * SIZE);
	}

	temp2 = (point*)malloc(sizeof(point) * w);
	memset(temp2, 0, sizeof(point) * w);
	for (i = 0; i < w; i++)
		temp2[i].flag = 0;

	d = (int)ceil(256 / (double)w);
	for (i = 0; i < w; i++)
	{
		j = 0;
		cnt = d * i;
		while (cnt >= 32)
		{
			cnt -= 32;
			j++;			//32bit가 넘으므로 32를 빼고 index를 1 증가
		}
		temp[i][j] = 1 << cnt;
		jacobianLtR(temp[i], P, p, &temp2[i]);		//temp[i]*P를 temp2[i]에 저장
	}

	/*	i = 2^4 + 2^2 + 2 + 1이라 하자
		precomputation[i] = 2^(4d)P + 2^(2d)P + 2^dP + P이다	*/
	for (i = 0; i < num; i++)
	{
		memset(&Q, NULL, sizeof(point));
		Q.flag = 1;
		cnt = 0;
		j = i;
		while (j)
		{
			if (j & 0x1)
				Q.flag = ECADD(&temp2[cnt], &Q, p, &Q);
			j >>= 1;
			cnt++;
		}
		memcpy(&precomputation[i], &Q, sizeof(point));
	}
}
//32bit word 8개로 구성된 256bit k를 w에 관한 comb에 맞게 2차원 배열 combK에 저장하는 함수
//combK는 w x d 배열이다.(d = 256/w) - 행 단위로 나눈다.
//combK[0] = {k(d-1), k(d-2), ... , k1, k0}
//k를 다 넣고 d의 배수로 떨어지지 않을 경우 MSB에 0padding을 필요한만큼 넣어준다.
void combPadding(IN word* k, IN int w, OUT char* combK[])
{
	int d, cnt1, cnt2, i;	//횟수를 증가시킬 cnt1, combK의 index를 위한 변수 cnt2
	word temp[SIZE] = { 0, }, temp2[SIZE] = { 0, };
	memcpy(temp, k, sizeof(word) * SIZE);

	d = ceil(256 / (double)w);
	cnt1 = 0, cnt2 = 0;
	for (cnt2 = 0; cnt2 < w; cnt2++)
	{
		for (cnt1 = 0; cnt1 < d; cnt1++)
		{
			combK[cnt2][cnt1] = temp[0] & 0x1;
			div2(temp);
		}
	}
}
void combSM(IN word* k, IN point* precomputation, IN word* p, IN char* combK[], IN int w, OUT point* kP)
{
	point Q;
	int i, j, num;
	memset(&Q, NULL, sizeof(point));
	Q.flag = 1;
	int d = ceil(256 / (double)w);

	for (i = d - 1; i >= 0; i--)
	{
		Q.flag = jacobianECDBL(&Q, p, &Q);
		//Q.flag = ECDBL(&Q, p, &Q);

		num = 0;
		for (j = 0; j < w; j++)
		{
			num += (combK[j][i] << j);
		}
		jacobianECADD(&Q, &precomputation[num], p, &Q);
		//ECADD(&Q, &precomputation[num], p, &Q);
	}
	//memcpy(kP, &Q, sizeof(point));
	jac2Aff(&Q, p, kP);
}
#endif
__int64 cpucycles(void)
{
	return __rdtsc();
}

int main()
{
	word opA[SIZE] = { 0x00000000, }, opB[SIZE] = { 0x00000000, }, opC[(2 * SIZE)] = { 0x00000000 , }, k[SIZE] = { 0, },
		p[SIZE] = { 0x00000000, }, modC[SIZE] = { 0x00000000, }, opC2[(2 * SIZE)] = { 0x00000000 , }, opInv[SIZE] = { 0, };

	byte buf[100];
	word** table;
	int w = 0, wNum = 0;
	point P, Q, R;
	point* precomputation;
	char** combK;
	//opA[0]이 최하위비트(제일 작은거), opA[7]이 최상위비트(제일 큰거)
	p[7] = 0xffffffff; p[6] = 0x00000001; p[5] = 0x00000000; p[4] = 0x00000000;
	p[3] = 0x00000000; p[2] = 0xffffffff; p[1] = 0xffffffff; p[0] = 0xffffffff;

	P.x[7] = 0x6b17d1f2; P.x[6] = 0xe12c4247; P.x[5] = 0xf8bce6e5; P.x[4] = 0x63a440f2;
	P.x[3] = 0x77037d81; P.x[2] = 0x2deb33a0; P.x[1] = 0xf4a13945; P.x[0] = 0xd898c296;

	P.y[7] = 0x4fe342e2; P.y[6] = 0xfe1a7f9b; P.y[5] = 0x8ee7eb4a; P.y[4] = 0x7c0f9e16;
	P.y[3] = 0x2bce3357; P.y[2] = 0x6b315ece; P.y[1] = 0xcbb64068; P.y[0] = 0x37bf51f5;

	P.flag = 0;

	FILE* fp1, * fp2, * fp3, * fp4;
	int i, j, num, d;
	ull cycles1 = 0, cycles2 = 0, cycles3 = 0, cycles4 = 0, cycles5 = 0, cycles6 = 0;

//CombMethod
#if 0
	fp1 = fopen("TV_Scalar.txt", "r");
	fp3 = fopen("ret_combSM.txt", "w");

	printf("원하는 w값을 입력하시오 : ");
	scanf("%d", &w);
	getchar();

	num = (int)pow(2, w);
	d = ceil(256 / (double)w);
	precomputation = (point*)malloc(sizeof(point) * num);
	precomputationComb(&P, w, num, p, precomputation);

	combK = (char**)malloc(sizeof(char*) * w);
	for (i = 0; i < w; i++)
	{
		combK[i] = (char*)malloc(sizeof(char) * d);
	}

	for (i = 0; i < 10000; i++)
	{
		readData(k, fp1);
		combPadding(k, w, combK);

		cycles1 = cpucycles();
		combSM(k, precomputation, p, combK, w, &Q);
		cycles2 = cpucycles();

		cycles3 += (cycles2 - cycles1);

		for (j = SIZE - 1; j >= 0; j--)
		{
			fprintf(fp3, "%08X", Q.x[j]);
		}
		fprintf(fp3, "\n");
		for (j = SIZE - 1; j >= 0; j--)
		{
			fprintf(fp3, "%08X", Q.y[j]);
		}
		fprintf(fp3, "\n\n");

		fgets(buf, sizeof(buf), fp1);
	}
	printf("%10lld\n", cycles3 / 10000);
#endif

//wNAFk기반 SM
#if 1
	word temp[SIZE] = { 1122334455,0,0,0,0,0,0,0 };	//k
	char wNAFk[257] = { 0, };						//wNAFk

	fp1 = fopen("TV_Scalar.txt", "r");
	fp3 = fopen("ret_wNAFSM.txt", "w");

	printf("원하는 w를 입력하시오 : ");
	scanf("%d", &w);

	table = (point**)malloc(sizeof(point*) * w);
	wNum = (int)pow(2, (double)w - 2);
	for (i = 0; i < wNum; i++)
		table[i] = (point*)malloc(sizeof(point) * SIZE);

	makeTable(&P, w, p, table);

	for (i = 0; i < 10000; i++)
	{
		readData(k, fp1);

		wNAF(k, w, wNAFk);
		cycles1 = cpucycles();
		wNAFSM(table, wNAFk, p, &Q);
		cycles2 = cpucycles();
		cycles3 += (cycles2 - cycles1);

		for (j = SIZE - 1; j >= 0; j--)
		{
			fprintf(fp3, "%08X", Q.x[j]);
		}
		fprintf(fp3, "\n");
		for (j = SIZE - 1; j >= 0; j--)
		{
			fprintf(fp3, "%08X", Q.y[j]);
		}
		fprintf(fp3, "\n\n");

		fgets(buf, sizeof(buf), fp1);
	}
	printf("%10lld\n", cycles3 / 10000);
#endif


//타원곡선연산 ECADD ECDBL jacobianECADD jacobianECDBL
#if 0
	fp1 = fopen("TV_Scalar.txt", "r");
	//fp2 = fopen("ret_SM.txt", "w");
	fp3 = fopen("ret_jacobian_SM.txt", "w");

	for (i = 0; i < 10000; i++)
	{
		readData(k, fp1);
		
		cycles1 = cpucycles();
		//LtR(k, &P, p, &Q);
		//RtL(k, &P, p, &Q);
		jacobianLtR(k, &P, p, &Q);
		cycles2 = cpucycles();
		cycles3 += (cycles2 - cycles1);

		//jacobian LtR
		for (j = SIZE - 1; j >= 0; j--)
		{
			fprintf(fp3, "%08X", Q.x[j]);
		}
		fprintf(fp3, "\n");
		for (j = SIZE - 1; j >= 0; j--)
		{
			fprintf(fp3, "%08X", Q.y[j]);
		}
		fprintf(fp3, "\n\n");

		//affine LtR RtL
		/*for (j = SIZE - 1; j >= 0; j--)
		{
			fprintf(fp2, "%08X", Q.x[j]);
		}
		fprintf(fp2, "\n");
		for (j = SIZE - 1; j >= 0; j--)
		{
			fprintf(fp2, "%08X", Q.y[j]);
		}
		fprintf(fp2, "\n\n");*/

		fgets(buf, sizeof(buf), fp1);
	}
	printf("%10lld\n", cycles3 / 10000);
#endif
//역원
#if 0
	fp1 = fopen("TV_opA_inv.txt", "r");
	fp3 = fopen("ret_PFINV.txt", "w");
	//fp4 = fopen("ret_SQR.txt", "w");

	for (i = 0; i < 10000; i++)
	{
		readData(opA, fp1);

		cycles1 = cpucycles();
		//fermatLT(opA, p, opInv);
		binaryInversion(opA, p, opInv);
		cycles2 = cpucycles();
		cycles3 += (cycles2 - cycles1);

		for (j = SIZE - 1; j >= 0; j--)
			fprintf(fp3, "%08X", opInv[j]);
		fprintf(fp3, "\n\n");

		fgets(buf, sizeof(buf), fp1);
}
	printf("%10lld\n", cycles3 / 10000);
#endif

//감산 
#if 0
	fp1 = fopen("TV_opA_PF.txt", "r");
	fp2 = fopen("TV_opB_PF.txt", "r");
	fp3 = fopen("ret_PFMUL.txt", "w");
	fp4 = fopen("ret_SQR.txt", "w");

	for (i = 0; i < 10000; i++)
	{
		readData(opA, fp1);
		readData(opB, fp2);

		multiplicationOS(opA, opB, opC);
		cycles1 = cpucycles();
		fastReduction(opC, p, modC);
		cycles2 = cpucycles();
		cycles3 += (cycles2 - cycles1);

		for (j = SIZE-1; j >= 0; j--)
			fprintf(fp3, "%08X", modC[j]);
		fprintf(fp3, "\n\n");

		fgets(buf, sizeof(buf), fp1);
		fgets(buf, sizeof(buf), fp2);
	}
	printf("%10lld\n", cycles3 / 10000);
#endif

//곱셈 OS, PS 연산 및 속도측정
#if 0
	fp1 = fopen("TV_opA_mul.txt", "r");
	fp2 = fopen("TV_opB_mul.txt", "r");
	fp3 = fopen("ret_MUL.txt", "w");
	fp4 = fopen("ret_SQR.txt", "w");

	for (i = 0; i < 10000; i++)
	{
		readData(opA, fp1);
		readData(opB, fp2);

		cycles1 = cpucycles();
		//multiplicationOS(opA, opA, opC);
		//multiplicationPS(opA, opB, opC);
		//mulOSdiv(opA, opB, opC);
		//mulPSdiv(opA, opB, opC);
		sqr(opA, opC2);
		cycles2 = cpucycles();
		cycles3 += (cycles2 - cycles1);

		for (j = ((2 * SIZE) - 1); j >= 0; j--)
			fprintf(fp3, "%08X", opC[j]);
		fprintf(fp3, "\n\n");

		/*for (j = ((2*SIZE)-1); j >= 0; j--)
			fprintf(fp4, "%08X", opC2[j]);
		fprintf(fp4, "\n\n");*/

		fgets(buf, sizeof(buf), fp1);
		fgets(buf, sizeof(buf), fp2);
	}
	printf("%10lld\n", cycles3 / 10000);
#endif
	
//덧셈 뺄셈 연산 및 속도측정
#if 0
	fp1 = fopen("TV_opA.txt", "r");
	fp2 = fopen("TV_opB.txt", "r");
	fp3 = fopen("ret_PFADD.txt", "w");
	fp4 = fopen("ret_PFSUB.txt", "w");

	for (j = 0; j < 10000; j++)
	{
		readData(opA, fp1);
		readData(opB, fp2);

		cycles1 = cpucycles();
		additionFp(opA, opB, p, opC);
		cycles2 = cpucycles();
		cycles3 += (cycles2 - cycles1);

		for (i = 7; i >= 0; i--)
			fprintf(fp3, "%08X", opC[i]);
		fprintf(fp3, "\n\n");

		cycles4 = cpucycles();
		subtractionFp(opA, opB, p, opC);
		cycles5 = cpucycles();
		cycles6 += (cycles5 - cycles4);

		for (i = 7; i >= 0; i--)
			fprintf(fp4, "%08X", opC[i]);
		fprintf(fp4, "\n\n");

		fgets(buf, sizeof(buf), fp1);
		fgets(buf, sizeof(buf), fp2);
	}

	printf("%10lld\n", cycles3 / 10000);
	printf("%10lld\n", cycles6 / 10000);
#endif

	return 0;
}
#endif
