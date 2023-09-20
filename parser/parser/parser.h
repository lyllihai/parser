#pragma once
//  PARSER.H
//  Header for Parser program.

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <string.h>
#pragma warning(disable:4996)
using namespace std;

//  declarations:
const double K = 1.38E-23;
const double Q = 1.60E-19;

double nodeValue[30] = { 0.0 }, jacMat[30][30] = { 0.0 }, result[30] = { 0.0 }, minDert[30] = { 0.0 }, initF[30] = { 0.0 }, preU = 0.0, stepSize = 0.0, stopTime = 0.0;
double initJac[30][30] = { 0.0 }, preX[30] = { 0.0 };
int Vsoure[10][4] = { 0 };    /*
						   Vsoure[x][0]; 表示V(x),是否链接两个非零节点
						   Vsoure[x][1]; 表示V(x),的F(x)是否已经输出一次
						   Vsoure[x][2]; 表示V(x),的第二个链接节点的nameNum
						   Vsoure[x][3]; 表示V(x),的jac矩阵是否已经输出一次

						  */


enum CompType {
	MOSFET, BJT, VSource, ISource, Inductor,
	Resistor, Diode, Capacitor
};
enum TranType { NMOS, PMOS, NPN, PNP };
enum Flag { UNSET, SET };
enum Boolean { FALSE, TRUE };
enum EquaType { Nodal, Modified };
int mCount = 0, bCount = 0, vSCount = 0, iSCount = 0, rCount = 0, iCount = 0, dCount = 0, cCount = 0, isTran = 0, stepNum = 0;
const int NameLength = 80, BufLength = 300, NA = -1;

class Component;
class ComponentHead;
class Node;
class NodeHead;
class Model;
class ModelHead;

struct Connectors {
	Flag flag;
	Node* node;
	int conNum;   //用户给节点的编号
};

struct Connections {
	Connections* next;
	Component* comp; // ~> pointer to component
	int conNum;
};


class Component {
public:
	Component(CompType typeIn, double valueIn, double tempIn, int con0In,
		int con1In, int con2In, int con3In, Model* modelIn, char* nameIn);
	~Component();
	CompType getType();
	int getNum();
	Component* getNext();
	void connect(int conNum, Node* nodeIn);
	double getVal();
	int getcompNum();
	void setNext(Component* nextIn);
	void setNum(int numIn);
	int getConVal(int conNum);
	Boolean isCon(int conNum);
	void print(int nodeNum, ofstream& outFile, int datum, int lastnode);
	void printMat(int nodeNum, int datum, int lastnode, double result[], int nameNum);
	void specialPrint(ofstream& outFile, int datum);
	void specialPrintMat(int datum, double result[]);
	void specialPrintJac(ofstream& outFile, int datum, Node* wrt/**/, int lastnode, EquaType eqType, Component* compPtr2, int* specPrintJacMNA /**/);
	void specialPrintJacMat(int datum, Node* wrt/**/, int lastnode, EquaType eqType, Component* compPtr2, int* specPrintJacMNA /**/, double jacMat[][30]);
	void printVal(ofstream& outFile);
	void printJac(int nodeNum, ofstream& outFile, int datum, int wrt, bool MNAflag);
	void printJacMat(int nodeNum, int datum, int wrt, bool MNAflag, double jacMat[][30], int fristIndex, int scendIndex);
	/*~> function to print the super node equation (for Nodal equation) */
	void printSuperNode(ofstream& outFile, int datum, int lastnode);
	void printSuperNodeMat(int datum, int lastnode, double result[]);
	/* */
	Node* getNode(int conNum);
	int getNodeNum(int conNum);
	char* getName();
private:
	Component* next;
	CompType type;
	Connectors con0, con1, con2, con3;
	int compNum;                                                            //同种类型组件的编号
	double value, temp;
	Model* model;
	char name[NameLength];
};

class Node {
public:
	Node(int num);
	~Node();
	int getNum();
	void setNameNum(int numIn);
	int getNameNum();
	int getCount();
	Connections* getConList();
	void connect(int conNumIn, Component* compIn);
	Node* getNext();
	void printNodal(ofstream& outFile, int datum, int lastnode);
	void printNodalMat(int datum, int lastnode, double result[]);
	void printMNA(ofstream& outFile, int datum, int lastnode);
	void printMNAMat(int datum, int lastnode, double result[]);
	void printJac(ofstream& outFile, int datum, Node* wrt, int lastnode, EquaType eqType);
	void printJacMat(int datum, Node* wrt, int lastnode, EquaType eqType, double jacMat[][30]);
	void printJacMNA(ofstream& outFile, int datum, Node* wrt, int lastnode);
	void printJacMNAMat(int datum, Node* wrt, int lastnode, double jacMat[][30]);
	void setNext(Node* nodeIn);
	/*~>function to print the super node equation (for Nodal equation)*/
	void printSuperNodal(ofstream& outFile, int datum, int lastnode);
	void printSuperNodalMat(int datum, int lastnode, double result[], int numIndex);
	/* */
private:
	Node* next;                                              //链接下一个节点，顺序和nodeList中的顺序保持一致
	int nodeNum, conCount;                                   //nodeNum 表示nodeList中节点的编号 从1开始，conCount表示该节点链接的边数
	Connections* conList;
	int nameNum;                                             //用户对节点的编号顺序
};

class NodeHead {
public:
	NodeHead();
	~NodeHead();
	Node* addNode();
	int getCount();
	Node* getNode(int nodeNum);
private:
	Node* nodeList;
	int nodeCount;
};

class CompHead {
public:
	CompHead();
	~CompHead();
	void addComp(Component* component);
	int getCount(CompType type);
	Component* getComp(int compNum);
private:
	Component* compList;
	int sCount, rCount, tCount, dCount, cCount;
};

class Model {
public:
	Model(char* nameIn, TranType typeIn, double isIn, double bfIn,
		double brIn, double tempIn);
	char* getName();
	TranType getType();
	double getIs();
	double getBf();
	double getBr();
	double getTemp();
	void setNext(Model* nextIn);
	Model* getNext();
private:
	char name[NameLength];
	double is, bf, br, temp;
	Model* next;
	TranType type;
};

class ModelHead {
public:
	ModelHead();
	void addModel(Model* modelIn);
	Model* getModel(char* nameIn);
private:
	Model* modelList;
};

//组件的构造函数
Component::Component(CompType typeIn, double valueIn = NA, double tempIn = NA,
	int con0In = NA, int con1In = NA, int con2In = NA, int con3In = NA,
	Model* modelIn = NULL, char* nameIn = NULL) {

	type = typeIn;
	con0.conNum = con0In;
	con1.conNum = con1In;
	con2.conNum = con2In;
	con3.conNum = con3In;
	con0.flag = UNSET;
	con1.flag = UNSET;
	con2.flag = UNSET;
	con3.flag = UNSET;
	value = valueIn;
	temp = tempIn;
	next = NULL;
	model = modelIn;
	strcpy(name, nameIn);

}
//组件的析构函数
Component::~Component() {};

void Component::connect(int conNum, Node* nodeIn) {
	//将元件的连接器与指定的节点进行连接。
	//通过传递连接器编号 conNum 和节点指针 nodeIn，函数将节点与相应的连接器进行连接。
	if (conNum == 0) {
		con0.node = nodeIn;
		con0.flag = SET;
	}
	if (conNum == 1) {
		con1.node = nodeIn;
		con1.flag = SET;
	}
	if (conNum == 2) {
		con2.node = nodeIn;
		con2.flag = SET;
	}
	if (conNum == 3) {
		con3.node = nodeIn;
		con3.flag = SET;
	}
}

CompType Component::getType() {
	return type;
}

int Component::getNum() {
	return compNum;
}

Component* Component::getNext() {
	return next;
}

double Component::getVal() {
	return value;
}

int Component::getcompNum() {
	return compNum;
}

void Component::setNext(Component* nextIn) {
	next = nextIn;
}

void Component::setNum(int numIn) {
	compNum = numIn;
}

int Component::getConVal(int conNum) {
	//返回指定连接器所连接的连接器编号
	int rtVal;
	if (conNum == 0)
		rtVal = con0.conNum;
	if (conNum == 1)
		rtVal = con1.conNum;
	if (conNum == 2)
		rtVal = con2.conNum;
	if (conNum == 3)
		rtVal = con3.conNum;
	return rtVal;
}

Boolean Component::isCon(int conNum) {
	//检查连接器是否连接
	Boolean rtVal;
	if (conNum == 0)
		rtVal = (con0.flag == SET) ? TRUE : FALSE;
	if (conNum == 1)
		rtVal = (con1.flag == SET) ? TRUE : FALSE;
	if (conNum == 2)
		rtVal = (con2.flag == SET) ? TRUE : FALSE;
	if (conNum == 3)
		rtVal = (con3.flag == SET) ? TRUE : FALSE;
	return rtVal;
}

void Component::print(int nodeNum, ofstream& outFile, int datum, int lastnode) {
	switch (type) {
	case MOSFET:
		outFile << name
			<< " connected: Drain->Node" << getNodeNum(0)
			<< " Source->Node" << getNodeNum(1)
			<< " Gate->Node" << getNodeNum(2)
			<< " Bulk->Node" << getNodeNum(3) << " ." << endl;


	case BJT:
		if ((con0.node->getNum() == nodeNum) && (model->getType() == NPN)) {
			outFile << " (" << name << "IS " << ")*(exp(-"
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")) -1)";
			outFile << " + ";
			outFile << "(-" << name << "IS "
				<< "/" << name << "BR "
				<< ")*(1+" << name << "BR "
				<< ")*(exp(-"
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")) -1) ";
		}

		if ((con2.node->getNum() == nodeNum) && (model->getType() == NPN)) {
			outFile << " (-" << name << "IS "
				<< "/" << name << "BF "
				<< ")*(1+" << name << "BF "
				<< ")*(exp(-"
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")) -1)";
			outFile << " - ";
			outFile << "(- " << name << "IS " << ")*(exp(-"
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")) -1) ";
		}

		if ((con1.node->getNum() == nodeNum) && (model->getType() == NPN)) {
			outFile << " (-" << name << "IS " << ")*(exp(-"
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")) -1)";
			outFile << " + ";
			outFile << "(" << name << "IS "
				<< "/" << name << "BR "
				<< ")*(1+" << name << "BR "
				<< ")*(exp(-"
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")) -1)";
			outFile << " + ";
			outFile << "(" << name << "IS "
				<< "/" << name << "BF "
				<< ")*(1+" << name << "BF "
				<< ")*(exp(-"
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")) -1)";
			outFile << " - ";
			outFile << "( " << name << "IS " << ")*(exp(-"
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")) -1) ";
		}

		if ((con0.node->getNum() == nodeNum) && (model->getType() == PNP)) {
			outFile << " (-" << name << "IS " << ")*(exp("
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")) -1)";
			outFile << " + ";
			outFile << "(" << name << "IS "
				<< "/" << name << "BR "
				<< ")*(1+" << name << "BR "
				<< ")*(exp("
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")) -1) ";
		}

		if ((con2.node->getNum() == nodeNum) && (model->getType() == PNP)) {
			outFile << " (" << name << "IS "
				<< "/" << name << "BF "
				<< ")*(1+" << name << "BF "
				<< ")*(exp("
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")) -1)";
			outFile << " - ";
			outFile << "( " << name << "IS " << ")*(exp("
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")) -1) ";
		}

		if ((con1.node->getNum() == nodeNum) && (model->getType() == PNP)) {
			outFile << " (" << name << "IS " << ")*(exp("
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")) -1)";
			outFile << " + ";
			outFile << "(-" << name << "IS "
				<< "/" << name << "BR "
				<< ")*(1+" << name << "BR "
				<< ")*(exp("
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")) -1)";
			outFile << " + ";
			outFile << "(-" << name << "IS "
				<< "/" << name << "BF "
				<< ")*(1+" << name << "BF "
				<< ")*(exp("
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")) -1)";
			outFile << " - ";
			outFile << "(-" << name << "IS " << ")*(exp("
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")) -1) ";
		}
		break;

	case VSource:
		if (con0.node->getNum() == nodeNum)
			outFile << " X(" << lastnode + compNum << ")";
		else if (con1.node->getNum() == nodeNum)
			outFile << " (-X(" << lastnode + compNum << ") )";
		break;

	case ISource:
		if (con0.node->getNum() == nodeNum)
			outFile << " (" << name << ") ";
		else if (con1.node->getNum() == nodeNum)
			outFile << " (-" << name << ") ";
		break;

	case Diode:
		if (con0.node->getNum() == nodeNum) {
			outFile << " (" << name << "IS " << ")*(exp("
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")) -1) ";
		}
		else if (con1.node->getNum() == nodeNum) {
			outFile << " (-" << name << "IS " << ")*(exp("
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")) -1) ";
		}
		break;

	case Resistor:
		if (con0.node->getNum() == nodeNum) {
			outFile << " (";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")/" << name << " ";
		}
		if (con1.node->getNum() == nodeNum) {
			outFile << " (";
			if (con1.node->getNameNum() != datum)
				outFile << "X(" << con1.node->getNameNum() << ')';
			if (con0.node->getNameNum() != datum)
				outFile << "-X(" << con0.node->getNameNum() << ')';
			outFile << ")/" << name << " ";
		}
		break;

	case Capacitor:
		outFile << " 0 ";
		break;

	case Inductor:
		if (con0.node->getNum() == nodeNum)
			outFile << " Il" << compNum << " ";
		else if (con1.node->getNum() == nodeNum)

			break;
	};
	return;
}

void Component::specialPrint(ofstream& outFile, int datum) {
	int eqnum;
	if (type == VSource) {
		if (con0.node->getNameNum() == datum) eqnum = con1.node->getNameNum();
		else eqnum = con0.node->getNameNum();
		outFile << endl << "F(" << eqnum << ") = ";
		outFile << " (";
		if (con0.node->getNameNum() != datum)
			outFile << "X(" << con0.node->getNameNum() << ')';
		if (con1.node->getNameNum() != datum)
			outFile << "-X(" << con1.node->getNameNum() << ')';
		outFile << ") -" << value << ';' << endl;
	}
	else if (type == Inductor) {
		outFile << endl << "F(I" << compNum << ") = ";
		outFile << " (";
		if (con0.node->getNameNum() != datum)
			outFile << "X(" << con0.node->getNameNum() << ')';
		if (con1.node->getNameNum() != datum)
			outFile << "-X(" << con1.node->getNameNum() << ')';
		outFile << ") " << ';' << endl;
	}
}

/// ~> POSSIBLE DEBUG: insert a condition with relation to the type of equation: MNA. i.e. if(typeEq = MNA) ...
/// ~> MODIFYING THE PROTOTYPE OF THE FUNCTION: including 'int lastnode' and 'EquaType eqType', Component* compPtr2

void Component::specialPrintJac(ofstream& outFile, int datum, Node* wrt /**/, int lastnode, EquaType eqType, Component* compPtr2, int* specPrintJacMNA /**/) {
	int eqnum;
	if (type == VSource) {
		if (con0.node->getNameNum() == datum)  eqnum = con1.node->getNameNum();
		else  eqnum = con0.node->getNameNum();
		outFile << endl << "JAC(" << eqnum << ", " << wrt->getNameNum() << ") = ";
		if (con0.node->getNameNum() == wrt->getNameNum())
			outFile << "1";
		else if (con1.node->getNameNum() == wrt->getNameNum())
			outFile << "-1";
		else
			outFile << "0";
		outFile << ';' << endl;
		//==================================================
		//~> Printing the jacobians of each voltage source with respect to each unknown current (MNEqs)
		if ((eqType == Modified) && (*specPrintJacMNA == 0)) {
			while (compPtr2 != NULL) {
				if (compPtr2->getType() == VSource) {

					outFile << endl << "JAC(" << eqnum << ", ";
					outFile << lastnode + compPtr2->getcompNum() << ") = 0; ";
					outFile << endl;

				}
				compPtr2 = compPtr2->getNext();
			}
			(*specPrintJacMNA) = 1;

		}
		//=====================================================

	}
	else if (type == Inductor) {
		outFile << endl << "JAC(I" << compNum << ", " << wrt->getNameNum() << ") = ";
		if (con0.node->getNameNum() == wrt->getNameNum())
			outFile << "1";
		else if (con1.node->getNameNum() == wrt->getNameNum())
			outFile << "-1";
		else
			outFile << "0";
		outFile << ';' << endl;
	}
}

void Component::printVal(ofstream& outFile) {
	switch (type) {
	case Resistor:
		outFile << name << " = " << value << ';' << endl;
		break;
	case VSource:
		outFile << name << " = " << value << ';' << endl;
		break;
	case Capacitor:
		outFile << name << " = " << value << ';' << endl;
		break;
	case BJT:
		outFile << name << "IS = " << model->getIs() << ';' << endl
			<< name << "BF = " << model->getBf() << ';' << endl
			<< name << "BR = " << model->getBr() << ';' << endl
			<< name << "N = " << (Q / (K * model->getTemp())) << ';' << endl;
		break;
	case Diode:
		outFile << name << "IS = " << value << ';' << endl
			<< name << "N = " << (Q / (K * temp)) << ';' << endl;
		break;
	};
}

//   	wrt = With Respect To
void Component::printJac(int nodeNum, ofstream& outFile, int datum, int wrt, bool MNAflag) {
	//~> MNAflag says if the type of analysis is MNA or NA
	switch (type) {
	case MOSFET:
		cout << "printJac for MOSFETs not implemented" << endl
			<< "PROGRAM ENDED ABNORMALLY!" << endl;
		exit(0);
		break;

	case BJT:
		if ((con0.node->getNum() == nodeNum) &&
			(model->getType() == NPN) &&
			(con0.node->getNameNum() == wrt)) {
			outFile << " (" << name << "IS"
				<< " * " << name << "N"
				<< " / " << name << "BR"
				<< ")*(1 + " << name << "BR"
				<< ")*(exp( - "
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}
		else if ((con0.node->getNum() == nodeNum) &&
			(model->getType() == NPN) &&
			(con2.node->getNameNum() == wrt)) {
			outFile << " ( - " << name << "IS"
				<< " * " << name << "N"
				<< ")*(exp( - "
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}
		else if ((con0.node->getNum() == nodeNum) &&
			(model->getType() == NPN) &&
			(con1.node->getNameNum() == wrt)) {
			outFile << " (" << name << "IS"
				<< " * " << name << "N"
				<< ")*(exp( - "
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
			outFile << " - ";
			outFile << "(" << name << "IS"
				<< " * " << name << "N"
				<< " / " << name << "BR"
				<< ")*(1 + " << name << "BR"
				<< ")*(exp( - "
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}

		else if ((con2.node->getNum() == nodeNum) &&
			(model->getType() == NPN) &&
			(con0.node->getNameNum() == wrt)) {
			outFile << " ( -  " << name << "IS"
				<< " * " << name << "N"
				<< ")*(exp(-"
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}
		else if ((con2.node->getNum() == nodeNum) &&
			(model->getType() == NPN) &&
			(con2.node->getNameNum() == wrt)) {
			outFile << " (" << name << "IS"
				<< " * " << name << "N"
				<< " / " << name << "BF"
				<< ")*(1 + " << name << "BF"
				<< ")*(exp(-"
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}
		else if ((con2.node->getNum() == nodeNum) &&
			(model->getType() == NPN) &&
			(con1.node->getNameNum() == wrt)) {
			outFile << " ( - " << name << "IS"
				<< " * " << name << "N"
				<< " / " << name << "BF"
				<< ")*(1 + " << name << "BF"
				<< ")*(exp(-"
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")))";
			outFile << " + ";
			outFile << "(" << name << "IS"
				<< " * " << name << "N"
				<< ")*(exp(-"
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}

		else if ((con1.node->getNum() == nodeNum) &&
			(model->getType() == NPN) &&
			(con0.node->getNameNum() == wrt)) {
			outFile << " ( - " << name << "IS"
				<< " * " << name << "N"
				<< " / " << name << "BR"
				<< ")*(1 + " << name << "BR"
				<< ")*(exp(-"
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << ")))";
			outFile << " + ";
			outFile << "(" << name << "IS"
				<< " * " << name << "N"
				<< ")*(exp(-"
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}
		else if ((con1.node->getNum() == nodeNum) &&
			(model->getType() == NPN) &&
			(con2.node->getNameNum() == wrt)) {
			outFile << " (" << name << "IS"
				<< " * " << name << "N"
				<< ")*(exp(-"
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << ")))";
			outFile << " - ";
			outFile << "(" << name << "IS"
				<< " * " << name << "N"
				<< " / " << name << "BF"
				<< ")*(1 + " << name << "BF"
				<< ")*(exp(-"
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}
		else if ((con1.node->getNum() == nodeNum) &&
			(model->getType() == NPN) &&
			(con1.node->getNameNum() == wrt)) {
			outFile << " ( - " << name << "IS"
				<< " * " << name << "N"
				<< ")*(exp(-"
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << ")))";
			outFile << " + ";
			outFile << "(" << name << "IS"
				<< " * " << name << "N"
				<< " / " << name << "BR "
				<< ")*(1 + " << name << "BR "
				<< ")*(exp(-"
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")))";
			outFile << " + ";
			outFile << "(" << name << "IS"
				<< " * " << name << "N"
				<< " / " << name << "BF "
				<< ")*(1 + " << name << "BF "
				<< ")*(exp(-"
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << ")))";
			outFile << " - ";
			outFile << "(" << name << "IS"
				<< " * " << name << "N"
				<< ")*(exp(-"
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}

		else if ((con0.node->getNum() == nodeNum) &&
			(model->getType() == PNP) &&
			(con0.node->getNameNum() == wrt)) {
			outFile << " (" << name << "IS"
				<< " * " << name << "N"
				<< " / " << name << "BR"
				<< ")*(1 + " << name << "BR"
				<< ")*(exp( "
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}
		else if ((con0.node->getNum() == nodeNum) &&
			(model->getType() == PNP) &&
			(con2.node->getNameNum() == wrt)) {
			outFile << " ( - " << name << "IS"
				<< " * " << name << "N"
				<< ")*(exp( "
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}
		else if ((con0.node->getNum() == nodeNum) &&
			(model->getType() == PNP) &&
			(con1.node->getNameNum() == wrt)) {
			outFile << " (" << name << "IS"
				<< " * " << name << "N"
				<< ")*(exp( "
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
			outFile << " - ";
			outFile << "(" << name << "IS"
				<< " * " << name << "N"
				<< " / " << name << "BR"
				<< ")*(1 + " << name << "BR"
				<< ")*(exp( "
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}

		else if ((con2.node->getNum() == nodeNum) &&
			(model->getType() == PNP) &&
			(con0.node->getNameNum() == wrt)) {
			outFile << " ( -  " << name << "IS"
				<< " * " << name << "N"
				<< ")*(exp("
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}
		else if ((con2.node->getNum() == nodeNum) &&
			(model->getType() == PNP) &&
			(con2.node->getNameNum() == wrt)) {
			outFile << " (" << name << "IS"
				<< " * " << name << "N"
				<< " / " << name << "BF"
				<< ")*(1 + " << name << "BF"
				<< ")*(exp("
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}
		else if ((con2.node->getNum() == nodeNum) &&
			(model->getType() == PNP) &&
			(con1.node->getNameNum() == wrt)) {
			outFile << " ( - " << name << "IS"
				<< " * " << name << "N"
				<< " / " << name << "BF"
				<< ")*(1 + " << name << "BF"
				<< ")*(exp("
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")))";
			outFile << " + ";
			outFile << "(" << name << "IS"
				<< " * " << name << "N"
				<< ")*(exp("
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}

		else if ((con1.node->getNum() == nodeNum) &&
			(model->getType() == PNP) &&
			(con0.node->getNameNum() == wrt)) {
			outFile << " ( - " << name << "IS"
				<< " * " << name << "N"
				<< " / " << name << "BR"
				<< ")*(1 + " << name << "BR"
				<< ")*(exp("
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << ")))";
			outFile << " + ";
			outFile << "(" << name << "IS"
				<< " * " << name << "N"
				<< ")*(exp("
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}
		else if ((con1.node->getNum() == nodeNum) &&
			(model->getType() == PNP) &&
			(con2.node->getNameNum() == wrt)) {
			outFile << " (" << name << "IS"
				<< " * " << name << "N"
				<< ")*(exp("
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << ")))";
			outFile << " - ";
			outFile << "(" << name << "IS"
				<< " * " << name << "N"
				<< " / " << name << "BF"
				<< ")*(1 + " << name << "BF"
				<< ")*(exp("
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}
		else if ((con1.node->getNum() == nodeNum) &&
			(model->getType() == PNP) &&
			(con1.node->getNameNum() == wrt)) {
			outFile << " ( - " << name << "IS"
				<< " * " << name << "N"
				<< ")*(exp("
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << ")))";
			outFile << " + ";
			outFile << "(" << name << "IS"
				<< " * " << name << "N"
				<< " / " << name << "BR "
				<< ")*(1 + " << name << "BR "
				<< ")*(exp("
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << "-X(" << con1.node->getNameNum() << ')';
			outFile << ")))";
			outFile << " + ";
			outFile << "(" << name << "IS"
				<< " * " << name << "N"
				<< " / " << name << "BF "
				<< ")*(1 + " << name << "BF "
				<< ")*(exp("
				<< name << "N*(";
			if (con2.node->getNameNum() != datum)
				outFile << "X(" << con2.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << ")))";
			outFile << " - ";
			outFile << "(" << name << "IS"
				<< " * " << name << "N"
				<< ")*(exp("
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}
		else
			outFile << " 0 ";
		break;

		//           NOTE:  Assuming DC voltages, the derivatives are all zero
	case VSource:
		if (((con0.node->getNum() == nodeNum) && (con0.node->getNameNum() == wrt) /*~>*/ && (MNAflag == FALSE)) ||
			((con1.node->getNum() == nodeNum) && (con1.node->getNameNum() == wrt) /*~>*/ && (MNAflag == FALSE)))
			outFile << " 1 ";
		else if (((con0.node->getNum() == nodeNum) && (con1.node->getNameNum() == wrt) /*~>*/ && (MNAflag == FALSE)) ||
			((con1.node->getNum() == nodeNum) && (con0.node->getNameNum() == wrt) /*~>*/ && (MNAflag == FALSE)))
			outFile << " (-1)";
		else
			outFile << " 0 ";
		break;

	case ISource:
		outFile << " 0 ";
		break;

	case Diode:
		if (((con0.node->getNum() == nodeNum) && (con0.node->getNameNum() == wrt)) ||
			((con1.node->getNum() == nodeNum) && (con1.node->getNameNum() == wrt))) {
			outFile << " (" << name << "IS"
				<< " * " << name << "N"
				<< ")*(exp("
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}
		else if (((con0.node->getNum() == nodeNum) && (con1.node->getNameNum() == wrt)) ||
			((con1.node->getNum() == nodeNum) && (con0.node->getNameNum() == wrt))) {
			outFile << " ( - " << name << "IS"
				<< " * " << name << "N"
				<< ")*(exp("
				<< name << "N*(";
			if (con0.node->getNameNum() != datum)
				outFile << "X(" << con0.node->getNameNum() << ')';
			if (con1.node->getNameNum() != datum)
				outFile << " - X(" << con1.node->getNameNum() << ')';
			outFile << "))) ";
		}
		else
			outFile << " 0 ";
		break;

	case Resistor:
		if (((con0.node->getNum() == nodeNum) && (con0.node->getNameNum() == wrt)) ||
			((con1.node->getNum() == nodeNum) && (con1.node->getNameNum() == wrt)))
			outFile << " 1/" << name << " ";
		else if (((con0.node->getNum() == nodeNum) && (con1.node->getNameNum() == wrt)) ||
			((con1.node->getNum() == nodeNum) && (con0.node->getNameNum() == wrt)))
			outFile << " (-1/" << name << ") ";
		else
			outFile << " 0 ";
		break;

	case Capacitor:
		outFile << " 0 ";
		break;

	case Inductor:
		cerr << "This section is not completed" << endl
			<< "PROGRAM ENDED ABNORMALLY!" << endl;
		exit(0);
		break;
	};
	return;
}

//节点法中，如果遇到独立电压源
void Component::printSuperNode(ofstream& outFile, int datum, int lastnode) {
	if (type == VSource) {
		if ((con0.node->getNameNum() != datum) && (con1.node->getNameNum() != datum))
		{ //~> if is a floating source
			outFile << endl << "F(" << con1.node->getNameNum() << ") = ";
			outFile << "(";
			con0.node->printSuperNodal(outFile, datum, lastnode);
			outFile << ") ";
			outFile << '+';
			outFile << " (";
			con1.node->printSuperNodal(outFile, datum, lastnode);
			outFile << ")";
			outFile << ';' << endl;
		}
	}
	return;
}

//获取组件所连接的节点的指针
Node* Component::getNode(int conNum) {
	switch (conNum) {
	case 0:
		return con0.node;
	case 1:
		return con1.node;
	case 2:
		return con2.node;
	case 3:
		return con3.node;
	}
}

//获取组件所连接的节点的编号
int Component::getNodeNum(int conNum) {
	switch (conNum) {
	case 0:
		return con0.node->getNum();
	case 1:
		return con1.node->getNum();
	case 2:
		return con2.node->getNum();
	case 3:
		return con3.node->getNum();
	}
}

char* Component::getName() {
	return name;
}

Node::Node(int Num) {
	next = NULL;
	nodeNum = Num;
	conCount = 0;
	conList = NULL;
	nameNum = 0;
}

Node::~Node() {};

int Node::getNum() {
	return nodeNum;
}

void Node::setNameNum(int numIn) {
	nameNum = numIn;
}

int Node::getNameNum() {
	return nameNum;
}

int Node::getCount() {
	return conCount;
}

Connections* Node::getConList() {
	return conList;
}

//将组件与节点相连接
void Node::connect(int conNumIn, Component* compIn) {
	Connections* conPtr;
	conCount++;
	if (conList == NULL) {
		conList = new Connections;
		conList->next = NULL;
		conList->conNum = conNumIn;//连接编号
		conList->comp = compIn;
	}
	else {
		conPtr = conList;
		while (conPtr->next != NULL)
			conPtr = conPtr->next;
		conPtr->next = new Connections;
		conPtr = conPtr->next;
		conPtr->next = NULL;
		conPtr->conNum = conNumIn;
		conPtr->comp = compIn;
	}
}

Node* Node::getNext() {
	return next;//返回下一节点
}

void Node::setNext(Node* nodeIn) {
	next = nodeIn;//设置下一节点
}

void Node::printNodal(ofstream& outFile, int datum, int lastnode) {
	Connections* conPtr;
	conPtr = conList;
	while (conPtr != NULL) { // ~> checks if some 'connection' is connected to any Voltage source. If there is no connections with Voltages source, print the equations.
		if (conPtr->comp->getType() == VSource)
			return;
		conPtr = conPtr->next;
	}
	outFile << endl << "F(" << nameNum << ") = ";
	conPtr = conList;
	while (conPtr->next != NULL) { //~> Using 'conPtr->next' not conPtr to check the end of the while loop
		conPtr->comp->print(nodeNum, outFile, datum, lastnode);
		outFile << '+';
		conPtr = conPtr->next;
	}
	conPtr->comp->print(nodeNum, outFile, datum, lastnode);
	outFile << ';' << endl;
	return;
}

void Node::printMNA(ofstream& outFile, int datum, int lastnode) {
	Connections* conPtr;
	int print;
	print = 0;//是否需要打印方程
	conPtr = conList;
	while (conPtr != NULL) {
		if (conPtr->comp->getType() == VSource) {  //~> seeking for a VSource
			if (Vsoure[conPtr->comp->getcompNum()][0] == 0) {
				print = 1;
				outFile << endl;
				outFile << "F(" << lastnode + conPtr->comp->getcompNum() << ")=";

			}
			else {
				if (Vsoure[conPtr->comp->getcompNum()][1] == 0) {
					print = 1;
					outFile << endl;
					outFile << "F(" << Vsoure[conPtr->comp->getcompNum()][2] << ")=";
					Vsoure[conPtr->comp->getcompNum()][1] = 1;
				}
				else {
					print = 1;
					outFile << endl;
					outFile << "F(" << lastnode + conPtr->comp->getcompNum() << ")=";
					Vsoure[conPtr->comp->getcompNum()][1] = 0;
				}
			}

		}
		conPtr = conPtr->next;
	}
	if (print) {
		conPtr = conList;
		while (conPtr->next != NULL) {
			conPtr->comp->print(nodeNum, outFile, datum, lastnode);
			outFile << "+";
			conPtr = conPtr->next;
		}
		conPtr->comp->print(nodeNum, outFile, datum, lastnode);
		outFile << ';' << endl;
	}
	return;
}

void Node::printJac(ofstream& outFile, int datum, Node* wrt, int lastnode, EquaType eqType) {
	Connections* conPtr;
	conPtr = conList;
	while (conPtr != NULL) {
		if (conPtr->comp->getType() == VSource)
			return;
		conPtr = conPtr->next;
	}
	outFile << endl << "JAC(" << nameNum << ", " << wrt->getNameNum() << ") = ";
	conPtr = conList;
	while (conPtr->next != NULL) {
		conPtr->comp->printJac(nodeNum, outFile, datum, wrt->getNameNum(), FALSE);
		outFile << "+";
		conPtr = conPtr->next;
	}
	conPtr->comp->printJac(nodeNum, outFile, datum, wrt->getNameNum(), FALSE);
	outFile << ';' << endl;
	if (eqType == Modified) {
		conPtr = wrt->conList;
		while (conPtr != NULL) {
			if (conPtr->comp->getType() == VSource) {
				if (Vsoure[conPtr->comp->getcompNum()][0] == 1) {
					if (Vsoure[conPtr->comp->getcompNum()][3] == 0) {
						outFile << endl << "JAC(" << nameNum << ", ";
						outFile << lastnode + conPtr->comp->getcompNum() << ") = 0; ";
						outFile << endl;

						Vsoure[conPtr->comp->getcompNum()][3] = 1;

					}
					else {
						Vsoure[conPtr->comp->getcompNum()][3] = 0;
					}
				}
				else {
					outFile << endl << "JAC(" << nameNum << ", ";
					outFile << lastnode + conPtr->comp->getcompNum() << ") = 0; ";
					outFile << endl;
				}


			}
			conPtr = conPtr->next;
		}
	}
	return;
}

void Node::printJacMNA(ofstream& outFile, int datum, Node* wrt, int lastnode) {
	Connections* conPtr;
	int print, val, tempVal;
	/* ~>*/
	bool MNAflagFALSE = FALSE; //Used to indicate if it is been printed JacMNA
	bool MNAflagTRUE = TRUE;
	/* <~*/
	print = 0;
	conPtr = conList;
	while (conPtr != NULL) {
		if (conPtr->comp->getType() == VSource) {

			if (Vsoure[conPtr->comp->getcompNum()][0] == 1) {
				if (conPtr->comp->getConVal(0) == nameNum) {
					print = 1;
					val = conPtr->comp->getcompNum();
					outFile << endl << "JAC(" << nameNum;
					outFile << ", " << wrt->getNameNum() << ") = ";
					tempVal = nameNum;
				}
				else {
					print = 1;
					val = conPtr->comp->getcompNum();
					outFile << endl << "JAC(" << lastnode + val;
					outFile << ", " << wrt->getNameNum() << ") = ";
					tempVal = lastnode + val;
				}
			}
			else {
				print = 1;
				val = conPtr->comp->getcompNum();
				outFile << endl << "JAC(" << lastnode + val;
				outFile << ", " << wrt->getNameNum() << ") = ";
				tempVal = lastnode + val;
			}




		}
		conPtr = conPtr->next;
	}
	if (print) {
		conPtr = conList;
		while (conPtr->next != NULL) {
			conPtr->comp->printJac(nodeNum, outFile, datum, wrt->getNameNum(), MNAflagTRUE);
			outFile << "+";
			conPtr = conPtr->next;
		}
		conPtr->comp->printJac(nodeNum, outFile, datum, wrt->getNameNum(), MNAflagTRUE);
		outFile << ';' << endl;
		conPtr = wrt->conList;
		while (conPtr != NULL) {
			if (conPtr->comp->getType() == VSource) {

				if (Vsoure[conPtr->comp->getcompNum()][0] == 1) {
					if (Vsoure[conPtr->comp->getcompNum()][3] == 0) {
						outFile << endl << "JAC(" << tempVal << ", ";
						outFile << lastnode + conPtr->comp->getcompNum() << ") = ";
						if (conPtr->comp->getcompNum() == val)  outFile << "1;";
						else outFile << "0;";
						outFile << endl;
						if (Vsoure[conPtr->comp->getcompNum()][0] == 1) {
							Vsoure[conPtr->comp->getcompNum()][3] = 1;
						}

					}
					else {
						Vsoure[conPtr->comp->getcompNum()][3] = 0;

					}
				}
				else {
					outFile << endl << "JAC(" << tempVal << ", ";
					outFile << lastnode + conPtr->comp->getcompNum() << ") = ";
					if (conPtr->comp->getcompNum() == val)  outFile << "1;";
					else outFile << "0;";
					outFile << endl;
				}



			}
			conPtr = conPtr->next;
		}
	}
	return;
}

void Node::printSuperNodal(ofstream& outFile, int datum, int lastnode) {
	Connections* conPtr;
	conPtr = conList;
	conPtr = conList;
	while (conPtr->next != NULL) { //~> Using 'conPtr->next' not conPtr to check the end of the while loop
		if (conPtr->comp->getType() != VSource) {
			conPtr->comp->print(nodeNum, outFile, datum, lastnode);
		}
		conPtr = conPtr->next;

	}
	if (conPtr->comp->getType() != VSource) {
		outFile << '+';
		conPtr->comp->print(nodeNum, outFile, datum, lastnode);
	}
	return;
}

NodeHead::NodeHead() {
	nodeList = NULL;
	nodeCount = 0;
}

NodeHead::~NodeHead() {};

Node* NodeHead::addNode() {
	Node* nodePtr;
	nodeCount++;
	if (nodeList == NULL) {
		nodeList = new Node(nodeCount);
		return nodeList;
	}
	else {
		nodePtr = nodeList;
		while (nodePtr->getNext() != NULL)
			nodePtr = nodePtr->getNext();
		nodePtr->setNext(new Node(nodeCount));
		return nodePtr->getNext();
	}
}


int NodeHead::getCount() {
	return nodeCount;
}

//获取指定节点
Node* NodeHead::getNode(int nodeNum) {
	Node* nodePtr;
	// need check that nodeNum does not exceed node count
	nodePtr = nodeList;
	for (int a = 0; a < nodeNum; a++)
		nodePtr = nodePtr->getNext();
	return nodePtr;
}

CompHead::CompHead() {
	compList = NULL;
	mCount = 0;
	bCount = 0;
	iCount = 0;
	rCount = 0;
	dCount = 0;
	cCount = 0;
	vSCount = 0;
	iSCount = 0;
}

CompHead::~CompHead() {};

void CompHead::addComp(Component* component) {
	Component* compPtr;
	switch (component->getType()) {
	case ISource:
		iSCount++;
		if (compList == NULL) {
			compList = component;
			component->setNum(iSCount);
		}
		else {
			compPtr = compList;
			while (compPtr->getNext() != NULL)
				compPtr = compPtr->getNext();
			compPtr->setNext(component);
			compPtr = compPtr->getNext();
			compPtr->setNum(iSCount);
		}
		break;
	case VSource:
		vSCount++;
		if (compList == NULL) {
			compList = component;
			component->setNum(vSCount);
		}
		else {
			compPtr = compList;
			while (compPtr->getNext() != NULL)
				compPtr = compPtr->getNext();
			compPtr->setNext(component);
			compPtr = compPtr->getNext();
			compPtr->setNum(vSCount);
		}
		break;
	case Resistor:
		rCount++;
		if (compList == NULL) {
			compList = component;
			component->setNum(rCount);
		}
		else {
			compPtr = compList;
			while (compPtr->getNext() != NULL)
				compPtr = compPtr->getNext();
			compPtr->setNext(component);
			compPtr = compPtr->getNext();
			compPtr->setNum(rCount);
		}
		break;
	case MOSFET:
		mCount++;
		if (compList == NULL) {
			compList = component;
			component->setNum(mCount);
		}
		else {
			compPtr = compList;
			while (compPtr->getNext() != NULL)
				compPtr = compPtr->getNext();
			compPtr->setNext(component);
			compPtr = compPtr->getNext();
			compPtr->setNum(mCount);
		}
		break;
	case BJT:
		bCount++;
		if (compList == NULL) {
			compList = component;
			component->setNum(bCount);
		}
		else {
			compPtr = compList;
			while (compPtr->getNext() != NULL)
				compPtr = compPtr->getNext();
			compPtr->setNext(component);
			compPtr = compPtr->getNext();
			compPtr->setNum(bCount);
		}
		break;
	case Diode:
		dCount++;
		if (compList == NULL) {
			compList = component;
			component->setNum(dCount);
		}
		else {
			compPtr = compList;
			while (compPtr->getNext() != NULL)
				compPtr = compPtr->getNext();
			compPtr->setNext(component);
			compPtr = compPtr->getNext();
			compPtr->setNum(dCount);
		}
		break;
	case Capacitor:
		cCount++;
		if (compList == NULL) {
			compList = component;
			component->setNum(cCount);
		}
		else {
			compPtr = compList;
			while (compPtr->getNext() != NULL)
				compPtr = compPtr->getNext();
			compPtr->setNext(component);
			compPtr = compPtr->getNext();
			compPtr->setNum(cCount);
		}
		break;
	case Inductor:
		iCount++;
		if (compList == NULL) {
			compList = component;
			component->setNum(iCount);
		}
		else {
			compPtr = compList;
			while (compPtr->getNext() != NULL)
				compPtr = compPtr->getNext();
			compPtr->setNext(component);
			compPtr = compPtr->getNext();
			compPtr->setNum(iCount);
		}
		break;
	}
}

int CompHead::getCount(CompType type) {
	switch (type) {
	case ISource:
		return iSCount;
	case VSource:
		return vSCount;
	case Resistor:
		return rCount;
	case Diode:
		return dCount;
	case MOSFET:
		return mCount;
	case BJT:
		return bCount;
	case Capacitor:
		return cCount;
	case Inductor:
		return iCount;
	}
}

Component* CompHead::getComp(int compNum) {
	Component* compPtr;
	compPtr = compList;
	for (int a = 0; a < compNum; a++)
		compPtr = compPtr->getNext();
	return compPtr;
}

Model::Model(char* nameIn, TranType typeIn, double isIn, double bfIn,
	double brIn, double tempIn) {
	strcpy(name, nameIn);
	type = typeIn;
	is = isIn;
	bf = bfIn;
	br = brIn;
	temp = tempIn;
	next = NULL;
}

TranType Model::getType() {
	return type;
}

char* Model::getName() {
	return name;
}

double Model::getIs() {
	return is;
}

double Model::getBf() {
	return bf;
}

double Model::getBr() {
	return br;
}

double Model::getTemp() {
	return temp;
}

void Model::setNext(Model* nextIn) {
	next = nextIn;
}

Model* Model::getNext() {
	return next;
}


ModelHead::ModelHead() {
	modelList = NULL;
}

void ModelHead::addModel(Model* modelIn) {
	Model* modelPtr;
	if (modelList == NULL) {
		modelList = modelIn;
	}
	else {
		modelPtr = modelList;
		while (modelPtr->getNext() != NULL) {
			modelPtr = modelPtr->getNext();
		}
		modelPtr->setNext(modelIn);
	}
}

Model* ModelHead::getModel(char* nameIn) {
	Model* modelPtr = modelList;
	while (strcmp(modelPtr->getName(), nameIn)) {
		modelPtr = modelPtr->getNext();
	}
	return modelPtr;
}






void Node::printNodalMat(int datum, int lastnode, double result[]) {
	Connections* conPtr;
	conPtr = conList;
	while (conPtr != NULL) { // ~> checks if some 'connection' is connected to any Voltage source. If there is no connections with Voltages source, print the equations.
		if (conPtr->comp->getType() == VSource)
			return;
		conPtr = conPtr->next;
	}
	conPtr = conList;
	while (conPtr->next != NULL) { //~> Using 'conPtr->next' not conPtr to check the end of the while loop
		conPtr->comp->printMat(nodeNum, datum, lastnode, result, nameNum);
		conPtr = conPtr->next;
	}
	conPtr->comp->printMat(nodeNum, datum, lastnode, result, nameNum);

	return;
}


void Component::printMat(int nodeNum, int datum, int lastnode, double result[], int nameNum) {
	switch (type) {


	case BJT:
		if ((con0.node->getNum() == nodeNum) && (model->getType() == NPN)) {
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + model->getIs() * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])) - 1);
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] + model->getIs() * (exp(-(Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]) - 1);
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + model->getIs() * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])) - 1);
			}


			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (-model->getIs() / model->getBr()) * (1 + model->getBr()) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])) - 1);
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] + (-model->getIs() / model->getBr()) * (1 + model->getBr()) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]) - 1);
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (-model->getIs() / model->getBr()) * (1 + model->getBr()) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])) - 1);
			}


		}

		if ((con2.node->getNum() == nodeNum) && (model->getType() == NPN)) {

			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (-model->getIs() / model->getBf()) * (1 + model->getBf()) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])) - 1);
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] + (-model->getIs() / model->getBf()) * (1 + model->getBf()) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]) - 1);
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (-model->getIs() / model->getBf()) * (1 + model->getBf()) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])) - 1);
			}



			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + model->getIs() * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])) - 1);
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] + model->getIs() * (exp(-(Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]) - 1);
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + model->getIs() * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])) - 1);
			}




		}

		if ((con1.node->getNum() == nodeNum) && (model->getType() == NPN)) {

			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] - model->getIs() * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])) - 1);
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] - model->getIs() * (exp(-(Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]) - 1);
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] - model->getIs() * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])) - 1);
			}


			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (model->getIs() / model->getBr()) * (1 + model->getBr()) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])) - 1);
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] + (model->getIs() / model->getBr()) * (1 + model->getBr()) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]) - 1);
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (model->getIs() / model->getBr()) * (1 + model->getBr()) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])) - 1);
			}



			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (model->getIs() / model->getBf()) * (1 + model->getBf()) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])) - 1);
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] + (model->getIs() / model->getBf()) * (1 + model->getBf()) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]) - 1);
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (model->getIs() / model->getBf()) * (1 + model->getBf()) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])) - 1);
			}




			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] - model->getIs() * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])) - 1);
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] - model->getIs() * (exp(-(Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]) - 1);
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] - model->getIs() * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])) - 1);
			}

		}

		if ((con0.node->getNum() == nodeNum) && (model->getType() == PNP)) {



			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] - model->getIs() * (exp((Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])) - 1);
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] - model->getIs() * (exp((Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]) - 1);
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] - model->getIs() * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])) - 1);
			}




			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (model->getIs() / model->getBr()) * (1 + model->getBr()) * (exp((Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])) - 1);
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] + (model->getIs() / model->getBr()) * (1 + model->getBr()) * (exp((Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]) - 1);
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (model->getIs() / model->getBr()) * (1 + model->getBr()) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])) - 1);
			}
		}

		if ((con2.node->getNum() == nodeNum) && (model->getType() == PNP)) {


			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (model->getIs() / model->getBf()) * (1 + model->getBf()) * (exp((Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])) - 1);
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] + (model->getIs() / model->getBf()) * (1 + model->getBf()) * (exp((Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]) - 1);
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (model->getIs() / model->getBf()) * (1 + model->getBf()) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])) - 1);
			}





			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] - model->getIs() * (exp((Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])) - 1);
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] - model->getIs() * (exp((Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]) - 1);
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] - model->getIs() * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])) - 1);
			}


		}

		if ((con1.node->getNum() == nodeNum) && (model->getType() == PNP)) {


			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + model->getIs() * (exp((Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])) - 1);
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] + model->getIs() * (exp((Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]) - 1);
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + model->getIs() * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])) - 1);
			}




			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (-model->getIs() / model->getBr()) * (1 + model->getBr()) * (exp((Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])) - 1);
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] + (-model->getIs() / model->getBr()) * (1 + model->getBr()) * (exp((Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]) - 1);
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (-model->getIs() / model->getBr()) * (1 + model->getBr()) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])) - 1);
			}




			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (-model->getIs() / model->getBf()) * (1 + model->getBf()) * (exp((Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])) - 1);
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] + (-model->getIs() / model->getBf()) * (1 + model->getBf()) * (exp((Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]) - 1);
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (-model->getIs() / model->getBf()) * (1 + model->getBf()) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])) - 1);
			}


			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + model->getIs() * (exp((Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])) - 1);
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] + model->getIs() * (exp((Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]) - 1);
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + model->getIs() * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])) - 1);
			}


		}
		break;

	case VSource:
		if (con0.node->getNum() == nodeNum)
			result[nameNum] = result[nameNum] + nodeValue[lastnode + compNum];
		else if (con1.node->getNum() == nodeNum)
			result[nameNum] = result[nameNum] - nodeValue[lastnode + compNum];
		break;

	case ISource:
		if (con0.node->getNum() == nodeNum)
			result[nameNum] = result[nameNum] + value;
		else if (con1.node->getNum() == nodeNum)
			result[nameNum] = result[nameNum] - value;
		break;

	case Diode:
		if (con0.node->getNum() == nodeNum) {


			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + model->getIs() * (exp((Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])) - 1);
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] + model->getIs() * (exp((Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]) - 1);
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + model->getIs() * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])) - 1);
			}
		}
		else if (con1.node->getNum() == nodeNum) {

			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] - model->getIs() * (exp((Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])) - 1);
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] - model->getIs() * (exp((Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]) - 1);
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] - model->getIs() * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])) - 1);
			}
		}
		break;

	case Resistor:
		if (con0.node->getNum() == nodeNum) {



			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()]) / value;
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] + (nodeValue[con0.node->getNameNum()]) / value;
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (-nodeValue[con1.node->getNameNum()]) / value;
			}
		}
		if (con1.node->getNum() == nodeNum) {


			if (con1.node->getNameNum() != datum && con0.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (nodeValue[con1.node->getNameNum()] - nodeValue[con0.node->getNameNum()]) / value;
			}
			if (con1.node->getNameNum() != datum && con0.node->getNameNum() == datum) {
				result[nameNum] = result[nameNum] + (nodeValue[con1.node->getNameNum()]) / value;
			}
			if (con1.node->getNameNum() == datum && con0.node->getNameNum() != datum) {
				result[nameNum] = result[nameNum] + (-nodeValue[con0.node->getNameNum()]) / value;
			}
		}
		break;

	case Capacitor:
		if (!isTran) {
			break;
		}
		else {
			result[nameNum] = result[nameNum] + value / stepSize * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()] - preX[con0.node->getNameNum()] + preX[con1.node->getNameNum()]);
			break;
		}

	case Inductor:
		if (con0.node->getNum() == nodeNum) {}

		else if (con1.node->getNum() == nodeNum)
		{
		}
		break;
	};
	return;
}


void Component::specialPrintMat(int datum, double result[]) {
	int eqnum;
	if (type == VSource) {
		if (con0.node->getNameNum() == datum) eqnum = con1.node->getNameNum();
		else eqnum = con0.node->getNameNum();


		if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
			result[eqnum] = result[eqnum] + (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()]) - value;
		}
		if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
			result[eqnum] = result[eqnum] + (nodeValue[con0.node->getNameNum()]) - value;
		}
		if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
			result[eqnum] = result[eqnum] + (-nodeValue[con1.node->getNameNum()]) - value;
		}
	}
	else if (type == Inductor) {
		//outFile << endl << "F(I" << compNum << ") = ";
		//outFile << " (";
		if (con0.node->getNameNum() != datum) {}
		//	outFile << "X(" << con0.node->getNameNum() << ')';
		if (con1.node->getNameNum() != datum) {}
		//outFile << "-X(" << con1.node->getNameNum() << ')';
	//outFile << ") " << ';' << endl;
	}
}


void Component::printSuperNodeMat(int datum, int lastnode, double result[]) {
	if (type == VSource) {
		if ((con0.node->getNameNum() != datum) && (con1.node->getNameNum() != datum))
		{

			int numIndex = con1.node->getNameNum();

			con0.node->printSuperNodalMat(datum, lastnode, result, numIndex);

			con1.node->printSuperNodalMat(datum, lastnode, result, numIndex);

		}
	}
	return;
}


void Node::printSuperNodalMat(int datum, int lastnode, double result[], int numIndex) {
	Connections* conPtr;
	conPtr = conList;
	conPtr = conList;
	while (conPtr->next != NULL) { //~> Using 'conPtr->next' not conPtr to check the end of the while loop
		if (conPtr->comp->getType() != VSource) {
			conPtr->comp->printMat(nodeNum, datum, lastnode, result, numIndex);
		}
		conPtr = conPtr->next;

	}
	if (conPtr->comp->getType() != VSource) {

		conPtr->comp->printMat(nodeNum, datum, lastnode, result, numIndex);
	}
	return;
}


void Node::printMNAMat(int datum, int lastnode, double result[]) {
	Connections* conPtr;
	int print, numIndex;
	print = 0;
	conPtr = conList;
	while (conPtr != NULL) {
		if (conPtr->comp->getType() == VSource) {  //~> seeking for a VSource

			if (Vsoure[conPtr->comp->getcompNum()][0] == 0) {
				print = 1;
				numIndex = lastnode + conPtr->comp->getcompNum();

			}
			else {
				if (Vsoure[conPtr->comp->getcompNum()][1] == 0) {
					print = 1;
					numIndex = Vsoure[conPtr->comp->getcompNum()][2];
					Vsoure[conPtr->comp->getcompNum()][1] = 1;

				}
				else {
					print = 1;
					numIndex = lastnode + conPtr->comp->getcompNum();
					Vsoure[conPtr->comp->getcompNum()][1] = 0;
				}
			}

		}


		conPtr = conPtr->next;
	}
	if (print) {
		conPtr = conList;
		while (conPtr->next != NULL) {
			conPtr->comp->printMat(nodeNum, datum, lastnode, result, numIndex);
			conPtr = conPtr->next;
		}
		conPtr->comp->printMat(nodeNum, datum, lastnode, result, numIndex);
	}
	return;
}


void Node::printJacMat(int datum, Node* wrt, int lastnode, EquaType eqType, double jacMat[][30]) {
	Connections* conPtr;
	int fristIndex = 0, scendIndex = 0;
	conPtr = conList;
	while (conPtr != NULL) {
		if (conPtr->comp->getType() == VSource)
			return;
		conPtr = conPtr->next;
	}

	fristIndex = nameNum; scendIndex = wrt->getNameNum();

	conPtr = conList;
	while (conPtr->next != NULL) {
		conPtr->comp->printJacMat(nodeNum, datum, wrt->getNameNum(), FALSE, jacMat, fristIndex, scendIndex);
		conPtr = conPtr->next;
	}
	conPtr->comp->printJacMat(nodeNum, datum, wrt->getNameNum(), FALSE, jacMat, fristIndex, scendIndex);
	if (eqType == Modified) {
		conPtr = wrt->conList;
		while (conPtr != NULL) {
			if (conPtr->comp->getType() == VSource) {
				if (Vsoure[conPtr->comp->getcompNum()][0] == 1) {
					if (Vsoure[conPtr->comp->getcompNum()][3] == 0) {
						jacMat[nameNum][lastnode + conPtr->comp->getcompNum()] = 0;

						Vsoure[conPtr->comp->getcompNum()][3] = 1;

					}
					else {
						Vsoure[conPtr->comp->getcompNum()][3] = 0;
					}
				}
				else {
					jacMat[nameNum][lastnode + conPtr->comp->getcompNum()] = 0;
				}

			}
			conPtr = conPtr->next;
		}
	}
	return;
}


void Component::printJacMat(int nodeNum, int datum, int wrt, bool MNAflag, double jacMat[][30], int fristIndex, int scendIndex) {
	//~> MNAflag says if the type of analysis is MNA or NA
	switch (type) {
	case MOSFET:
		cout << "printJac for MOSFETs not implemented" << endl
			<< "PROGRAM ENDED ABNORMALLY!" << endl;
		exit(0);
		break;

	case BJT:
		if ((con0.node->getNum() == nodeNum) &&
			(model->getType() == NPN) &&
			(con0.node->getNameNum() == wrt)) {


			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]));
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}
		}
		else if ((con0.node->getNum() == nodeNum) &&
			(model->getType() == NPN) &&
			(con2.node->getNameNum() == wrt)) {

			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]));
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}

		}
		else if ((con0.node->getNum() == nodeNum) &&
			(model->getType() == NPN) &&
			(con1.node->getNameNum() == wrt)) {

			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]));
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}



			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] - (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] - (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]));
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] - (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}
		}

		else if ((con2.node->getNum() == nodeNum) &&
			(model->getType() == NPN) &&
			(con0.node->getNameNum() == wrt)) {


			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]));
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}
		}
		else if ((con2.node->getNum() == nodeNum) &&
			(model->getType() == NPN) &&
			(con2.node->getNameNum() == wrt)) {


			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]));
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}
		}
		else if ((con2.node->getNum() == nodeNum) &&
			(model->getType() == NPN) &&
			(con1.node->getNameNum() == wrt)) {

			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]));
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}

			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]));
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}
		}

		else if ((con1.node->getNum() == nodeNum) &&
			(model->getType() == NPN) &&
			(con0.node->getNameNum() == wrt)) {


			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] - (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] - (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]));
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] - (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}


			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]));
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}
		}
		else if ((con1.node->getNum() == nodeNum) &&
			(model->getType() == NPN) &&
			(con2.node->getNameNum() == wrt)) {

			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]));
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}

			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]));
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}
		}
		else if ((con1.node->getNum() == nodeNum) &&
			(model->getType() == NPN) &&
			(con1.node->getNameNum() == wrt)) {

			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]));
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}


			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]));
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}


			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]));
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}


			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]));
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp(-(Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}
		}

		else if ((con0.node->getNum() == nodeNum) &&
			(model->getType() == PNP) &&
			(con0.node->getNameNum() == wrt)) {


			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp((Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp((Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]));
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}
		}
		else if ((con0.node->getNum() == nodeNum) &&
			(model->getType() == PNP) &&
			(con2.node->getNameNum() == wrt)) {

			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]));
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}
		}
		else if ((con0.node->getNum() == nodeNum) &&
			(model->getType() == PNP) &&
			(con1.node->getNameNum() == wrt)) {


			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]));
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}


			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] - (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp((Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] - (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp((Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]));
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] - (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}
		}

		else if ((con2.node->getNum() == nodeNum) &&
			(model->getType() == PNP) &&
			(con0.node->getNameNum() == wrt)) {

			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]));
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}
		}
		else if ((con2.node->getNum() == nodeNum) &&
			(model->getType() == PNP) &&
			(con2.node->getNameNum() == wrt)) {


			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp((Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp((Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]));
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}
		}
		else if ((con2.node->getNum() == nodeNum) &&
			(model->getType() == PNP) &&
			(con1.node->getNameNum() == wrt)) {

			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp((Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp((Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]));
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}


			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]));
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}
		}

		else if ((con1.node->getNum() == nodeNum) &&
			(model->getType() == PNP) &&
			(con0.node->getNameNum() == wrt)) {

			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] - (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp((Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] - (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp((Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]));
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] - (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}


			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]));
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}
		}
		else if ((con1.node->getNum() == nodeNum) &&
			(model->getType() == PNP) &&
			(con2.node->getNameNum() == wrt)) {

			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]));
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}


			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp((Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp((Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]));
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}
		}
		else if ((con1.node->getNum() == nodeNum) &&
			(model->getType() == PNP) &&
			(con1.node->getNameNum() == wrt)) {

			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]));
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}


			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp((Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp((Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]));
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBr()) * (1 + model->getBr()) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}


			if (con2.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp((Q / (K * model->getTemp())) * (nodeValue[con2.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con2.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp((Q / (K * model->getTemp())) * nodeValue[con2.node->getNameNum()]));
			}
			if (con2.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp())) / model->getBf()) * (1 + model->getBf()) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}


			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]));
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}
		}
		else
			jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + 0;
		break;

		//           NOTE:  Assuming DC voltages, the derivatives are all zero
	case VSource:
		if (((con0.node->getNum() == nodeNum) && (con0.node->getNameNum() == wrt) /*~>*/ && (MNAflag == FALSE)) ||
			((con1.node->getNum() == nodeNum) && (con1.node->getNameNum() == wrt) /*~>*/ && (MNAflag == FALSE)))
			jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + 1;
		else if (((con0.node->getNum() == nodeNum) && (con1.node->getNameNum() == wrt) /*~>*/ && (MNAflag == FALSE)) ||
			((con1.node->getNum() == nodeNum) && (con0.node->getNameNum() == wrt) /*~>*/ && (MNAflag == FALSE)))
			jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] - 1;
		else
			jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + 0;
		break;

	case ISource:
		jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + 0;
		break;

	case Diode:
		if (((con0.node->getNum() == nodeNum) && (con0.node->getNameNum() == wrt)) ||
			((con1.node->getNum() == nodeNum) && (con1.node->getNameNum() == wrt))) {

			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]));
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}
		}
		else if (((con0.node->getNum() == nodeNum) && (con1.node->getNameNum() == wrt)) ||
			((con1.node->getNum() == nodeNum) && (con0.node->getNameNum() == wrt))) {

			if (con0.node->getNameNum() != datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (nodeValue[con0.node->getNameNum()] - nodeValue[con1.node->getNameNum()])));
			}
			if (con0.node->getNameNum() != datum && con1.node->getNameNum() == datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * nodeValue[con0.node->getNameNum()]));
			}
			if (con0.node->getNameNum() == datum && con1.node->getNameNum() != datum) {
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + (-model->getIs() * (Q / (K * model->getTemp()))) * (exp((Q / (K * model->getTemp())) * (-nodeValue[con1.node->getNameNum()])));
			}
		}
		else
			jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + 0;
		break;

	case Resistor:
		if (((con0.node->getNum() == nodeNum) && (con0.node->getNameNum() == wrt)) ||
			((con1.node->getNum() == nodeNum) && (con1.node->getNameNum() == wrt)))
			jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + 1 / value;

		else if (((con0.node->getNum() == nodeNum) && (con1.node->getNameNum() == wrt)) ||
			((con1.node->getNum() == nodeNum) && (con0.node->getNameNum() == wrt)))
			jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] - 1 / value;
		else
			jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + 0;
		break;

	case Capacitor:
		if (!isTran) {
			jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + 0;
			break;
		}
		else {
			if (((con0.node->getNum() == nodeNum) && (con0.node->getNameNum() == wrt)) ||
				((con1.node->getNum() == nodeNum) && (con1.node->getNameNum() == wrt)))
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + value / stepSize;

			else if (((con0.node->getNum() == nodeNum) && (con1.node->getNameNum() == wrt)) ||
				((con1.node->getNum() == nodeNum) && (con0.node->getNameNum() == wrt)))
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] - value / stepSize;
			else
				jacMat[fristIndex][scendIndex] = jacMat[fristIndex][scendIndex] + 0;
			break;
		}


	case Inductor:
		cerr << "This section is not completed" << endl
			<< "PROGRAM ENDED ABNORMALLY!" << endl;
		exit(0);
		break;
	};
	return;
}


void Component::specialPrintJacMat(int datum, Node* wrt /**/, int lastnode, EquaType eqType, Component* compPtr2, int* specPrintJacMNA /**/, double jacMat[][30]) {
	int eqnum;
	if (type == VSource) {
		if (con0.node->getNameNum() == datum)  eqnum = con1.node->getNameNum();
		else  eqnum = con0.node->getNameNum();

		if (con0.node->getNameNum() == wrt->getNameNum())
			jacMat[eqnum][wrt->getNameNum()] = jacMat[eqnum][wrt->getNameNum()] + 1;

		else if (con1.node->getNameNum() == wrt->getNameNum())
			jacMat[eqnum][wrt->getNameNum()] = jacMat[eqnum][wrt->getNameNum()] - 1;
		else
			jacMat[eqnum][wrt->getNameNum()] = jacMat[eqnum][wrt->getNameNum()] + 0;
		//==================================================
		//~> Printing the jacobians of each voltage source with respect to each unknown current (MNEqs)
		if ((eqType == Modified) && (*specPrintJacMNA == 0)) {
			while (compPtr2 != NULL) {
				if (compPtr2->getType() == VSource) {
					jacMat[eqnum][lastnode + compPtr2->getcompNum()] = 0;
				}
				compPtr2 = compPtr2->getNext();
			}
			(*specPrintJacMNA) = 1;

		}
		//=====================================================

	}
	else if (type == Inductor) {
		//outFile << endl << "JAC(I" << compNum << ", " << wrt->getNameNum() << ") = ";
		if (con0.node->getNameNum() == wrt->getNameNum()) {}
		//outFile << "1";
		else if (con1.node->getNameNum() == wrt->getNameNum()) {}
		//outFile << "-1";
		else {}
		//outFile << "0";
	//outFile << ';' << endl;
	}
}


void Node::printJacMNAMat(int datum, Node* wrt, int lastnode, double jacMat[][30]) {
	Connections* conPtr;
	int print, val, fristIndex, scendIndex, tempVal;
	/* ~>*/
	bool MNAflagFALSE = FALSE; //Used to indicate if it is been printed JacMNA
	bool MNAflagTRUE = TRUE;
	/* <~*/
	print = 0;
	conPtr = conList;
	while (conPtr != NULL) {
		if (conPtr->comp->getType() == VSource) {

			if (Vsoure[conPtr->comp->getcompNum()][0] == 1) {
				if (conPtr->comp->getConVal(0) == nameNum) {
					print = 1;
					val = conPtr->comp->getcompNum();
					fristIndex = nameNum;
					scendIndex = wrt->getNameNum();
					tempVal = nameNum;
				}
				else {

					print = 1;

					val = conPtr->comp->getcompNum();
					fristIndex = lastnode + val;
					scendIndex = wrt->getNameNum();
					tempVal = lastnode + val;
				}
			}
			else {
				print = 1;
				val = conPtr->comp->getcompNum();
				fristIndex = lastnode + val;
				scendIndex = wrt->getNameNum();
				tempVal = lastnode + val;
			}


		}
		conPtr = conPtr->next;
	}
	if (print) {
		conPtr = conList;
		while (conPtr->next != NULL) {
			conPtr->comp->printJacMat(nodeNum, datum, wrt->getNameNum(), MNAflagTRUE, jacMat, fristIndex, scendIndex);
			conPtr = conPtr->next;
		}
		conPtr->comp->printJacMat(nodeNum, datum, wrt->getNameNum(), MNAflagTRUE, jacMat, fristIndex, scendIndex);
		conPtr = wrt->conList;
		while (conPtr != NULL) {
			if (conPtr->comp->getType() == VSource) {
				if (Vsoure[conPtr->comp->getcompNum()][0] == 1) {
					if (Vsoure[conPtr->comp->getcompNum()][3] == 0) {
						if (conPtr->comp->getcompNum() == val)
							jacMat[tempVal][lastnode + conPtr->comp->getcompNum()] = 1;
						else
							jacMat[tempVal][lastnode + conPtr->comp->getcompNum()] = 0;
						if (Vsoure[conPtr->comp->getcompNum()][0] == 1) {
							Vsoure[conPtr->comp->getcompNum()][3] = 1;


						}

					}
					else {
						Vsoure[conPtr->comp->getcompNum()][3] = 0;

					}
				}
				else {
					if (conPtr->comp->getcompNum() == val)
						jacMat[tempVal][lastnode + conPtr->comp->getcompNum()] = 1;
					else
						jacMat[tempVal][lastnode + conPtr->comp->getcompNum()] = 0;
				}

			}
			conPtr = conPtr->next;
		}
	}
	return;
}




























