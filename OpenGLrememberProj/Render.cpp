

#include "Render.h"

#include <Windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>
#include <thread>





inline double f_1(double p0, double p1, double p2, double p3,  double t)
{
	return p0 * pow((1 - t), 3) + 3 * t * pow((1 - t), 2) * p1 + 3 * t * t*(1 - t) * p2 + p3 * t * t * t; //формула Безье
	
}

void curve_Bezier(double P1[], double P2[], double P3[], double P4[], double delta_time) {
	
	glBegin(GL_LINE_STRIP);
	for (double t = 0; t <= 1; t += 0.001) {
		double P[3];
		P[0] = f_1(P1[0], P2[0], P3[0], P4[0], t);
		P[1] = f_1(P1[1], P2[1], P3[1], P4[1], t);
		P[2] = f_1(P1[2], P2[2], P3[2], P4[2], t);
		glVertex3dv(P);
	}
	glEnd();
	//Наши векторы для 1 кривой	
	glColor3b(0, 0, 1);
	glLineWidth(1);
	
	glBegin(GL_LINE_STRIP);
	glVertex3dv(P1);
	glVertex3dv(P2);
	glVertex3dv(P3);
	glVertex3dv(P4);
	glEnd();

}

double f_2(double p1, double p4, double r1, double r4, double t)
{
	return p1 * (2 * t * t * t - 3 * t * t + 1) + p4 * (3 * t * t - 2 * t * t * t) + r1 * (t * t * t - 2 * t * t + t) + r4 * (t * t * t - t * t); //Формула Эрмита 
}

void curve_Hermite(double P1[], double P2[], double P3[], double P4[]) {
	glPopMatrix();
	glPushMatrix(); 
	glLineWidth(3);
	glColor3d(1, 0, 0);
	glBegin(GL_LINE_STRIP);
	for (double t = 0; t <= 1; t += 0.001) {
		double P[3];
		P[0] = f_2(P1[0], P2[0], P3[0], P4[0], t);
		P[1] = f_2(P1[1], P2[1], P3[1], P4[1], t);
		P[2] = f_2(P1[2], P2[2], P3[2], P4[2], t);
		glVertex3dv(P);
	}
	glEnd();
	glLineWidth(1.5);
	glColor3d(1, 0, 1);
	glBegin(GL_LINE_STRIP);
	glVertex3dv(P1);
	glVertex3dv(P3);
	glEnd();
	glLineWidth(1.5);
	glBegin(GL_LINE_STRIP);
	glVertex3dv(P2);
	glVertex3dv(P4);
	glEnd();

	glPointSize(6.5);
	glColor3d(0, 0.2, 0.4);
	glBegin(GL_POINTS);
	glVertex3dv(P1);
	glVertex3dv(P2);
	glEnd();
	glPopMatrix();

}

double factorial(double num) {
	double res=1;
	for (int i = 1; i <= num; i++) {
		res *= i;
	}
	return res;
}

double Bernstein_Polynomial(double n, double i, double u) {
	return ((factorial(n) / (factorial(i) * factorial(n - i))) * pow(u, i) * pow((1 - u), (n - i)));
}

double* Bezier_Surface_calc(double P[4][4][3], double u, double v) {

	double P123[3];
	P123[0] = 0;
	P123[1] = 0;
	P123[2] = 0;
	double P2[3] = { P[1][3][0],P[1][3][1],P[1][3][2] };
	for (int i = 0; i <= 3; i++) {
		for (int j = 0; j <= 3; j++) {
			P123[0] += Bernstein_Polynomial(3, i, u) * Bernstein_Polynomial(3, j, v) * P[i][j][0];
			P123[1] += Bernstein_Polynomial(3, i, u) * Bernstein_Polynomial(3, j, v) * P[i][j][1];
			P123[2] += Bernstein_Polynomial(3, i, u) * Bernstein_Polynomial(3, j, v) * P[i][j][2];
		}
	}
	return P123;
}

void Bezier_Surface() {

	double* Trace = new double[3];
	glPushMatrix();
	glTranslated(9, 4, 0);
	double P[4][4][3] = {
		{{0,9,1}, {3,11,0}, {6,9,0}, {10,9,1}},
		{{0,7,0}, {3,6,0}, {6,6,4}, {9,6,0}},
		{{4,3,0}, {3,4,3}, {7,3,1}, {9,4,0}},
		{{0,0,1}, {3,0,5}, {9,0,0}, {12,0,1}},
	};
	int i = 0;
	int j = 0;
	glPointSize(6);
	glColor3d(0, 1, 0);
	glBegin(GL_POINTS); // заданные точки, соединенные между собой
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			glVertex3dv(P[i][j]);
		}
		glVertex3dv(P[i][j]);
	}
	glEnd();

	i = 0;
	j = 0;
	glColor3d(0, 0, 1);

	while (i < 4) {
		glBegin(GL_LINE_STRIP); //соединение, между заданными точками
		while (j < 4) {
			glVertex3dv(P[i][j]);
			j += 1;
		}
		j = 0;
		i += 1;
		glEnd();
	}
	i = 0;
	j = 0;


	while (i < 4) {
		glBegin(GL_LINE_STRIP); //соединение, между заданными точками
		while (j < 4) {
			glVertex3dv(P[j][i]);
			j += 1;
		}
		j = 0;
		i += 1;
		glEnd();
	}

	glTranslated(0, 12, 0);
	glPointSize(6.5);
	glColor3d(1, 1, 0.2);

	glBegin(GL_POINTS);  // поверхность Безъе точечная

	i = 0;
	j = 0;
	for (double u = 0; u <= 1; u += 0.1) {
		j = 0;
		for (double v = 0; v <= 1; v += 0.1) {
			Trace = Bezier_Surface_calc(P, u, v);

			glVertex3dv(Trace);

			j += 1;
		}
		i += 1;
	}
	glEnd();
	glPopMatrix();

	glTranslated(9, -6, 0);
	glColor3d(1, 1, 0);


	double u = 0;
	double v = 0;
	glColor3d(1, 0.5, 0.5);

	while (u <= 1) {
		glBegin(GL_LINE_STRIP); // поверхность Безье из линий
		while (v <= 1) {
			Trace = Bezier_Surface_calc(P, u, v);

			glVertex3dv(Trace);
			v += 0.1;
		}
		v = 0;
		u += 0.1;
		glEnd();
	}

	u = 0;
	v = 0;
	while (u <= 1) {
		glBegin(GL_LINE_STRIP);
		while (v <= 1) {
			Trace = Bezier_Surface_calc(P, v, u);

			glVertex3dv(Trace);
			v += 0.1;
		}
		v = 0;
		u += 0.1;
		glEnd();
	}
}

double  PointsDistance(double P1[3], double P2[3]) {
	return sqrt(pow(P2[0] - P1[0], 2) + pow(P2[1] - P1[0], 2) + pow(P2[2] - P1[2], 2));
}


double FindAngle(double V1[3], double V2[3]) {
	double c = V1[0] * V2[0] + V1[1] * V2[1] + V1[2] * V2[2];
	double d = sqrt(pow(V1[0], 2) + pow(V1[1], 2) + pow(V1[2], 2)) * sqrt(pow(V2[0], 2) + pow(V2[1], 2) + pow(V2[2], 2));
	return c / d;
}

double t_max = 0;

double direction = 1;

void PlaneCurse(double P1[], double P4[], double R1[], double R4[], double t_max) {

	double Vec[3] = { f_2(P1[0], P4[0], R1[0], R4[0], t_max) ,f_2(P1[1], P4[1], R1[1], R4[1], t_max),f_2(P1[2], P4[2], R1[2], R4[2], t_max) }; //постороение точек кривой
	double Vec1[3] = { f_2(P1[0], P4[0], R1[0], R4[0], t_max + 0.01 * direction) ,f_2(P1[1], P4[1], R1[1], R4[1], t_max + 0.01 * direction),f_2(P1[2], P4[2], R1[2], R4[2], t_max + 0.01 * direction) }; //постороение точек кривой +1	
	double dirV[3] = { Vec1[0] - Vec[0], Vec1[1] - Vec[1], Vec1[2] - Vec[2] }; //вектор между двумя точками	
	double nulXYZ[3] = { 0, 0, 0 };

	double length = PointsDistance(nulXYZ, dirV);
	double normal_length = (1 / length); //нормализация вектора
	double x = dirV[0] * normal_length;
	double y = dirV[1] * normal_length;
	double z = dirV[2] * normal_length;
	double XYZ[3] = { x,y,z }; //нормализованный вектор
	double orig[3] = { 1,0,0 }; //
	double rotatedXY[3] = { XYZ[0],XYZ[1],0 }; //кручение по ХУ
	
	length = PointsDistance(nulXYZ, rotatedXY);
	normal_length = (1 / length); //нормализация направления кручения
	x = rotatedXY[0] * normal_length;
	y = rotatedXY[1] * normal_length;
	z = rotatedXY[2] * normal_length;
	double rotatedX[3] = { x,y,z };
	double cosXY = FindAngle(orig, rotatedX); //нахождение угла между векторами
	double vecpr[3] = { orig[1] * rotatedX[2] - orig[2] * rotatedX[1], orig[2] * rotatedX[0] - orig[0] * rotatedX[2], orig[0] * rotatedX[1] - orig[1] * rotatedX[0] };
	double sinSign = vecpr[2] / abs(vecpr[2]);
	double angleXY = acos(cosXY) * 180 / acos(-1) * sinSign;
	double origZ[3] = { 0, 0, 1 };
	double cosXZ = FindAngle(origZ, XYZ);
	double angleXZ = (acos(cosXZ) * 180.0 / acos(-1)) - 90;
	

	glPushMatrix();
	glTranslated(Vec[0], Vec[1], Vec[2]);
	glRotated(angleXY, 0, 0, 1); //ху
	glRotated(angleXZ, 0, 1, 0); //z

	glRotated(-90, 0, 0, 1);


	glBegin(GL_TRIANGLES);
	glVertex3f(-1, 0.0, 0);
	glVertex3f(0, 1.5, 0);
	glVertex3f(1, 0, 0);
	glEnd();

	glBegin(GL_TRIANGLES);
	glVertex3f(-0.5, 0.75, 0);
	glVertex3f(-1, 0.75, 0);
	glVertex3f(-0.75, 0.25, 0);
	glEnd();

	glBegin(GL_TRIANGLES);
	glVertex3f(0.5, 0.75, 0);
	glVertex3f(1, 0.75, 0);
	glVertex3f(0.75, 0.25, 0);
	glEnd();
	
	glColor3d(1, 0.5, 0);
	glBegin(GL_QUADS);
	glVertex3f(-1, 0.75, 0);
	glVertex3f(-1, 1.0, 0);
	glVertex3f(-0.5, 1.0, 0);
	glVertex3f(-0.5, 0.75, 0);
	glEnd();

	glColor3d(1, 0.5, 0);
	glBegin(GL_QUADS);
	glVertex3f(1, 0.75, 0);
	glVertex3f(1, 1.0, 0);
	glVertex3f(0.5, 1.0, 0);
	glVertex3f(0.5, 0.75, 0);
	glEnd();

	glPopMatrix();
}



void Render(double delta_time)
{    
	t_max += delta_time / 5 * direction; //t_max становится = 1 за 5 секунд
	if (t_max > 1) direction = -1;
	if (t_max < 0) direction = 1;


	//Безье
	//Наши точки для 1 кривой	
	double P0_1[] = { 0,0,0 };
	double P1_1[] = { 1,-2,0 }; 
	double P2_1[] = { -4,6,7 };
	double P3_1[] = { 10,10,0 };

	//Наши точки для 2 кривой
	double P0_2[] = { 0,0,0 };
	double P1_2[] = { -6,4,0 };
	double P2_2[] = { -5,6,7 };
	double P3_2[] = { 10,-10,4 };
	//Наши точки для Эрмита
	glColor3b(90, 0, 0);
	//curve_Bezier(P0_1, P1_1, P2_1, P3_1, delta_time);
	glColor3b(0, 90, 0);
	curve_Bezier(P0_2, P1_2, P2_2, P3_2, delta_time);


	
	double P1_3[] = { 0,0,0 }; 
	double P2_3[] = { 5,4,8 };
	double P3_3[] = { 2,7,-6 };
	double P4_3[] = { 2, 8,-9 };

	double P1_4[] = { 0,0,0};
	double P2_4[] = { 3,1,10};
	double R3_4[] = { -3,5,-4};
	double R4_4[] = { 0,3,-5};
	curve_Hermite(P1_3, P2_3, P3_3, P4_3);
	PlaneCurse(P1_3, P2_3, P3_3, P4_3, t_max);

	Bezier_Surface();
}   

