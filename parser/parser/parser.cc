#include "parser.h"
using namespace std;
using namespace Eigen;
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <string>

double stripString(char* stringIn);
void printComponents(Component* compPtr);
void printNodes(Node* nodePtr, int compFlag);
char* strComponentType(Component* compPtr);
char* ComponentTypeName(Component* compPtr);  //obtain component type name
int portNum(Component* comPtr, Node* nodePtr); //obtain port number
bool isAccurate(double result[], int num, double accurateValue);
void NR_Iterations(double jacMat[][30], double result[], double minDert[], int number, int& count, double accurateValue, int datum, int lastnode, bool Homotopy = false, double t = 0);
void get_circuit_equation(int datum, int lastnode);
void get_Jacobian_matrix(int datum, int lastnode);
void Homotopy_NR(double jacMat[][30], double result[], double minDert[], int number, int& count, double accurateValue, int datum, int lastnode, double t);




NodeHead nodeList;
CompHead compList;

int main(int argc, char* argv[]) {
    ifstream inFile;
    ofstream outFile;
    ofstream outfile;       //add another 'outfile' to finish relevant work
    ModelHead modelList;

    // Buffers used in parsing:
    char inName[NameLength], outName[NameLength], buf[BufLength], myOutName[NameLength],
        buf1[BufLength], buf2[BufLength], buf3[BufLength], nameBuf[NameLength],
        * bufPtr, * charPtr1, * charPtr2;
    int intBuf1, intBuf2, intBuf3, intBuf4, datum = NA, eqNum = NA, specPrintJacMNA = 0;
    double douBuf1, douBuf2, douBuf3, douBuf4;
    CompType typeBuf;
    Component* compPtr, * compPtr1, * compPtr2;
    Node* nodePtr, * nodePtr1, * nodePtr2;
    Model* modelPtr;
    TranType TtypeBuf;
    EquaType eqType = Modified;

    strcpy(inName, "NOTHING");
    strcpy(outName, "NOTHING");
    strcpy(myOutName, "NOTHING");

    // process equation types:
    if (eqNum == NA) {
        while ((eqNum != 1) && (eqNum != 2)) {
            cout << "Available Equations Types Are:" << endl
                << " <1>  Nodal" << endl
                << " <2>  Modified Nodal" << endl
                << "Please enter your choice <1, 2>:" << endl;
            cin >> buf;
            eqNum = atoi(buf);
        }
        if (eqNum == 1)
            eqType = Nodal;
        else if (eqNum == 2)
            eqType = Modified;
    }

    // process input file name:
    if (!strcmp(inName, "NOTHING")) {
        cerr << "Please enter the input Spice Netlist: <QUIT to exit>" << endl;
        cin >> inName;
        if (!strcmp(inName, "QUIT")) {
            cerr << "Program Exited Abnormally!" << endl;
            exit(0);
        }
    }
    inFile.open(inName, ios::in);
    while (!inFile) {
        cerr << inName << " is an invalid input file." << endl
            << "Please enter the input Spice Netlist: <QUIT to exit>" << endl;
        cin >> inName;
        if (!strcmp(inName, "QUIT")) {
            cerr << "Program Exited Abnormally!" << endl;
            exit(0);
        }
        inFile.open(inName, ios::in);
    }

    // process output file
    if (!strcmp(outName, "NOTHING")) {
        strcpy(outName, inName);
        strtok(outName, ".");
        strcat(outName, ".Pout");
    }
    outFile.open(outName, ios::out);
    cout << endl;


    //***************读取网表内容*************
    inFile.getline(buf, BufLength);       // first line of netlist is discarded
    inFile.getline(buf, BufLength);


    while (inFile.good()) {
        if ((buf == NULL) || (*buf == '\0')) {
            inFile.getline(buf, BufLength);
            continue;
        }
        strcpy(buf1, buf);
        if (!strcmp(strtok(buf1, " "), ".model")) {
            strcpy(buf2, strtok(NULL, " "));
            charPtr1 = strtok(NULL, " ");
            if (!strcmp(charPtr1, "PNP"))
                TtypeBuf = PNP;
            else if (!strcmp(charPtr1, "NPN"))
                TtypeBuf = NPN;
            else if (!strcmp(charPtr1, "NMOS"))
                TtypeBuf = NMOS;
            else if (!strcmp(charPtr1, "PMOS"))
                TtypeBuf = PMOS;

            charPtr1 = strtok(NULL, " ");
            while (charPtr1 != NULL) {
                strcpy(buf3, "");
                if ((charPtr1[0] == 'I') && (charPtr1[1] == 'S') && (charPtr1[2] == '=')) {
                    douBuf1 = stripString(charPtr1);
                }
                if ((charPtr1[0] == 'B') && (charPtr1[1] == 'F') && (charPtr1[2] == '=')) {
                    douBuf2 = stripString(charPtr1);
                }
                if ((charPtr1[0] == 'B') && (charPtr1[1] == 'R') && (charPtr1[2] == '=')) {
                    douBuf3 = stripString(charPtr1);
                }
                if ((charPtr1[0] == 'T') && (charPtr1[1] == 'E') && (charPtr1[2] == '=')) {
                    douBuf4 = stripString(charPtr1);
                }
                charPtr1 = strtok(NULL, " ");
            }
            modelPtr = new Model(buf2, TtypeBuf, douBuf1, douBuf2, douBuf3, douBuf4);
            modelList.addModel(modelPtr);
        }
        inFile.getline(buf, BufLength);
    }
    inFile.close();
    inFile.open(inName, ios::in);


    //.Tran
    inFile.getline(buf, BufLength);       // first line of netlist is discarded
    inFile.getline(buf, BufLength);


    while (inFile.good()) {
        if ((buf == NULL) || (*buf == '\0')) {
            inFile.getline(buf, BufLength);
            continue;
        }
        strcpy(buf1, buf);
        if (!strcmp(strtok(buf1, " "), ".tran")) {

            charPtr1 = strtok(NULL, " ");
            stepSize = atof(charPtr1);

            charPtr1 = strtok(NULL, " ");
            stopTime = atof(charPtr1);
        }
        inFile.getline(buf, BufLength);
    }
    inFile.close();
    inFile.open(inName, ios::in);


    char model_str[9];
    //  starting of parsing by creating linked list of components
    inFile.getline(buf, BufLength);       // first line of netlist is discarded
    inFile.getline(buf, BufLength);
    while (inFile.good()) {
        if ((buf == NULL) || (*buf == '\0')) {
            inFile.getline(buf, BufLength);
            continue;
        }

        if (isalpha(*buf)) {

            //  EDIT THIS SECTION IF NEW COMPONENTS ARE ADDED!!!
            //  we could do some rearranging in this section to catch each type in order.
            switch (*buf) {
            case 'v':
            case 'V':
            {
                int vsnum = vSCount;
                typeBuf = VSource;
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                douBuf1 = atof(strtok(NULL, " "));
                compPtr = new Component(typeBuf, douBuf1, NA, intBuf1, intBuf2, NA, NA, NULL, nameBuf);
                compList.addComp(compPtr);
                if (intBuf1 != 0 && intBuf2 != 0) {
                    vsnum++;
                    Vsoure[vsnum][0] = 1;
                    Vsoure[vsnum][2] = intBuf2;
                }
                break;
            }

            case 'i':
            case 'I':
                cout << "I" << endl;
                typeBuf = ISource;
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                douBuf1 = atof(strtok(NULL, " "));
                compPtr = new Component(typeBuf, douBuf1, NA, intBuf1, intBuf2, NA, NA, NULL, nameBuf);
                compList.addComp(compPtr);
                break;
            case 'q':
            case 'Q':
                typeBuf = BJT;
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                intBuf3 = atoi(strtok(NULL, " "));
                compPtr = new Component(typeBuf, NA, NA, intBuf1, intBuf2, intBuf3, NA,
                    modelList.getModel(strtok(NULL, " ")), nameBuf);
                compList.addComp(compPtr);
                break;
            case 'm':
            case 'M':
                typeBuf = MOSFET;
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                intBuf3 = atoi(strtok(NULL, " "));
                intBuf4 = atoi(strtok(NULL, " "));
                compPtr = new Component(typeBuf, NA, NA, intBuf1, intBuf2, intBuf3, intBuf4,
                    modelList.getModel(strtok(NULL, " ")), nameBuf);
                compList.addComp(compPtr);
                break;
            case 'r':
            case 'R':
                typeBuf = Resistor;
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                douBuf1 = atof(strtok(NULL, " "));
                compPtr = new Component(typeBuf, douBuf1, NA, intBuf1, intBuf2, NA, NA, NULL, nameBuf);
                compList.addComp(compPtr);
                break;
            case 'd':
            case 'D':
                typeBuf = Diode;
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                charPtr1 = strtok(NULL, " ");
                while (charPtr1 != NULL) {
                    if ((charPtr1[0] == 'I') && (charPtr1[1] == 'S') && (charPtr1[2] == '=')) {
                        douBuf1 = stripString(charPtr1);
                    }
                    if ((charPtr1[0] == 'T') && (charPtr1[1] == 'E') && (charPtr1[4] == '=')) {
                        douBuf2 = stripString(charPtr1);
                    }
                    charPtr1 = strtok(NULL, " ");
                }
                compPtr = new Component(typeBuf, douBuf1, douBuf2, intBuf1, intBuf2, NA, NA, NULL, nameBuf);
                compList.addComp(compPtr);
                break;
            case 'c':
            case 'C':
                typeBuf = Capacitor;
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                douBuf1 = atof(strtok(NULL, " "));
                compPtr = new Component(typeBuf, douBuf1, NA, intBuf1, intBuf2, NA, NA, NULL, nameBuf);
                compList.addComp(compPtr);
                break;
            case 'l':
            case 'L':
                typeBuf = Inductor;
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                douBuf1 = atof(strtok(NULL, " "));
                compPtr = new Component(typeBuf, douBuf1, NA, intBuf1, intBuf2, NA, NA, NULL, nameBuf);
                compList.addComp(compPtr);
                break;
            };
        }
        inFile.getline(buf, BufLength);
    }


    //  Now the components are created and it is time to set up the list of nodes.
    //  we should actually use second connector of first Source as the first Node (Datum)
    compPtr1 = compList.getComp(0);
    while (compPtr1 != NULL) {
        for (int b = 0; b < 3; b++) { /* ~> J. Erik Melo note: A component can have until 4 connectors. But here just 3 are been considered. It should change the condition to 'b <= 3' or 'b < 4'?*/
            if ((!compPtr1->isCon(b)) && (compPtr1->getConVal(b) != NA)) { //~> verify if the connector 'b' is not set && if the name of the node to which this same connector 'b' is connected is a valid name as found in the circuit file. That is, if the name is not NA, that is, if this connector was named in the instantiation of the component.
                intBuf1 = compPtr1->getConVal(b); // ~> getting the connector number as in the netlist file
                nodePtr1 = nodeList.addNode();
                nodePtr1->setNameNum(intBuf1);  // ~> naming the node as in the netlist file
                compPtr1->connect(b, nodePtr1); // ~> connecting the 'connector' of component to the node
                nodePtr1->connect(b, compPtr1); // ~> connecting the 'connection' of the node to the component

                // now search and connect all other appropriate connectors to this node.
                // error checking should be added to prevent duplicated, or skipped connectors.
                compPtr2 = compPtr1->getNext();
                while (compPtr2 != NULL) {
                    for (int c = 0; c < 3; c++) { //~> verifying which one of the others connectors (of components) are connected to the node above
                        if (compPtr2->getConVal(c) == intBuf1) { //~> if next component in the list of components has a connector with the same name (conNum) of the connector above, connect it to the same node.
                            compPtr2->connect(c, nodePtr1);
                            nodePtr1->connect(c, compPtr2);
                            break;                                    //~> As a component can only have one connector with the same name (connected in the same node), don't search the others and go out of the 'for' loop
                        }
                    }
                    compPtr2 = compPtr2->getNext();
                }
            }
        }
        compPtr1 = compPtr1->getNext();
    }


    //  At this point, we are done creating a representation of the circuit in memory
    //  now, we need to call each node to create and output its nodal equation.
    //  Each node will call the components attached for the individual contributions to the
    //  nodal equation.

      // verify that input datum is valid
    Boolean check = FALSE;
    if (datum != NA) {
        nodePtr = nodeList.getNode(0);
        while (nodePtr != NULL) {
            if (nodePtr->getNameNum() == datum)
                check = TRUE;
            nodePtr = nodePtr->getNext();
        }
        if (check == FALSE) {
            cerr << "Datum value invalid!" << endl
                << "PROGRAM EXITED ABNORMALLY!" << endl;
            exit(0);
        }
    }

    // Loop to find lastnode
    //找到节点列表中的最后一个节点，并获取其名称编号。
    nodePtr = nodeList.getNode(0); //~> getting the pointer to the first node, pointed by 'headNode'
    int lastnode = nodePtr->getNameNum();
    while (nodePtr != NULL) {
        lastnode = (nodePtr->getNameNum() > lastnode) ? nodePtr->getNameNum() : lastnode;
        nodePtr = nodePtr->getNext();
    }

    //  Loop to find the datum
    if (datum == NA) {
        nodePtr = nodeList.getNode(0);
        nodePtr1 = nodePtr->getNext();
        while (nodePtr1 != NULL) {
            if (nodePtr1->getCount() > nodePtr->getCount())
                nodePtr = nodePtr1;
            nodePtr1 = nodePtr1->getNext();
        }
        //datum = nodePtr->getNameNum();
        datum = 0; //此处做出了修改，暂时不明白该处理有啥作用
    }

    //=================================
    //~> Checking the component list
    //~> Comment this part to omit
    compPtr = compList.getComp(0);
    printComponents(compPtr);

    nodePtr = nodeList.getNode(0);
    printNodes(nodePtr, 1);

    //<~
    //==================================

    // output circuit information
    outFile << "%Parser V1.0" << endl;
    outFile << "%Input Spice Deck:  " << inName << endl;
    outFile << "%Equation Type:     ";
    if (eqType == Nodal)
        outFile << "NODAL" << endl;
    else if (eqType == Modified)
        outFile << "MODIFIED NODAL" << endl;
    outFile << "%Datum Node:        " << datum << endl;


    // create value table
    outFile << endl
        << "%*****************************************************************************" << endl;
    outFile << "%                      Component Values:" << endl;
    compPtr = compList.getComp(0);
    while (compPtr != NULL) {
        compPtr->printVal(outFile);
        compPtr = compPtr->getNext();
    }
    outFile << endl
        << "%*****************************************************************************" << endl;

    // go down the nodal list and have components announce themselves
    outFile << endl << "%                      Circuit Equations: " << endl;
    nodePtr = nodeList.getNode(0);
    while (nodePtr != NULL) {
        if (nodePtr->getNameNum() != datum) {
            nodePtr->printNodal(outFile, datum, lastnode);
        }
        nodePtr = nodePtr->getNext();
    }

    //go down the component list and give equations for all sources
    compPtr = compList.getComp(0);
    while (compPtr != NULL) {
        compPtr->specialPrint(outFile, datum);
        compPtr = compPtr->getNext();
    }

    //~> go down the component list and give supernode equations for all float sources (Nodal Analysis)
    if (eqType != Modified) {
        compPtr = compList.getComp(0);
        while (compPtr != NULL) {
            compPtr->printSuperNode(outFile, datum, lastnode);
            compPtr = compPtr->getNext();
        }
    }


    // go down the node list and give additional MNA equations
    if (eqType == Modified) {
        nodePtr = nodeList.getNode(0);
        while (nodePtr != NULL) {
            if (nodePtr->getNameNum() != datum)
                nodePtr->printMNA(outFile, datum, lastnode);
            nodePtr = nodePtr->getNext();
        }
    }


    // print jacobians
    outFile << endl
        << "%*****************************************************************************" << endl;
    outFile << endl << "%                      Jacobians: " << endl;
    nodePtr1 = nodeList.getNode(0);
    while (nodePtr1 != NULL) {   //~> this loop handles the nodes not connected to a Vsource and those ones that are not the 'datum' node
        if (nodePtr1->getNameNum() != datum) {
            nodePtr2 = nodeList.getNode(0);
            while (nodePtr2 != NULL) {
                if (nodePtr2->getNameNum() != datum) {
                    nodePtr1->printJac(outFile, datum, nodePtr2, lastnode, eqType);
                }
                nodePtr2 = nodePtr2->getNext();
            }
        }
        nodePtr1 = nodePtr1->getNext();
    }

    // go down the component list and give equations for all sources
    compPtr = compList.getComp(0);
    while (compPtr != NULL) {
        nodePtr2 = nodeList.getNode(0);
        compPtr2 = compList.getComp(0);
        while (nodePtr2 != NULL) {
            if (nodePtr2->getNameNum() != datum) {
                compPtr->specialPrintJac(outFile, datum, nodePtr2, lastnode, eqType, compPtr2, &specPrintJacMNA); // ~> specPrintJacMNA is used to verify if the jacobians w.r.t. the Modified equations was already printed to print only once.
            }
            nodePtr2 = nodePtr2->getNext();
        }
        specPrintJacMNA = 0;
        compPtr = compPtr->getNext();
    }


    // print the Jacobians for the additional MNA equations
    if (eqType == Modified) {
        nodePtr1 = nodeList.getNode(0);
        while (nodePtr1 != NULL) {
            if (nodePtr1->getNameNum() != datum) {
                nodePtr2 = nodeList.getNode(0);
                while (nodePtr2 != NULL) {
                    if (nodePtr2->getNameNum() != datum)

                        nodePtr1->printJacMNA(outFile, datum, nodePtr2, lastnode);
                    nodePtr2 = nodePtr2->getNext();
                }
            }
            nodePtr1 = nodePtr1->getNext();
        }
    }


    cout << endl;


    // %***************************************************************************************************

    if (!strcmp(myOutName, "NOTHING")) {
        strcpy(myOutName, inName);
        strtok(myOutName, ".");
        strcat(myOutName, "out.txt");
    }
    outfile.open(myOutName, ios::out);


    nodePtr = nodeList.getNode(0);

    outfile << "datum = " << datum << "\t\t" << "lastnode = " << lastnode << endl;
    Connections* conPtr;
    while (nodePtr != NULL) {

        outfile << "节点 " << nodePtr->getNameNum() << "\t\t" << "所连器件数为：" << nodePtr->getCount() << endl;
        conPtr = nodePtr->getConList();
        while (conPtr->next != NULL) {
            outfile << "\t\t" << "编号： " << conPtr->comp->getcompNum() << "\t\t" << "类型： " << ComponentTypeName(conPtr->comp) << "\t\t" << "链接端口：" << portNum(conPtr->comp, nodePtr) << "\t\t";
            switch (conPtr->comp->getType()) {
            case VSource:
                outfile << "名称：" << "VCC" << endl; break;
            default:
                outfile << "名称：" << strComponentType(conPtr->comp) << conPtr->comp->getcompNum() << endl; break;
            }
            outfile << "\t\t" << "value:" << conPtr->comp->getVal() << endl;

            conPtr = conPtr->next;
        }

        outfile << "\t\t" << "编号： " << conPtr->comp->getcompNum() << "\t\t" << "类型： " << ComponentTypeName(conPtr->comp) << "\t\t" << "链接端口：" << portNum(conPtr->comp, nodePtr) << "\t\t";
        switch (conPtr->comp->getType()) {
        case VSource:
            outfile << "名称：" << "VCC" << endl; break;
        default:
            outfile << "名称：" << strComponentType(conPtr->comp) << conPtr->comp->getcompNum() << endl;
            break;
        }
        outfile << "\t\t" << "value:" << conPtr->comp->getVal() << endl;

        nodePtr = nodePtr->getNext();
    }

    outfile.close();

    //%****************************************************************************************************
    int choose = 0;

    cout << "*************************************" << endl;

    cout << "   1、N-R" << endl;
    cout << "   2、Homotopy" << endl;
    cout << "   3、Transient（瞬态、伪瞬态）" << endl;

    cout << "*************************************" << endl;

    cin >> choose;
    
    if (choose == 1) {
        //以下实现矩阵的运算

        int number = 0;
        int count = 1;
        double accurateValue;

        cout << "Please enter the initial data number：" << endl;
        cin >> number;
        cout << "Please enter the initial data value:" << endl;
        for (int i = 0; i < number; i++) {
            cin >> nodeValue[i + 1];
        }
        cout << "please input required accuracy:" << endl;
        cin >> accurateValue;
        
        NR_Iterations(jacMat, result, minDert, number, count, accurateValue, datum, lastnode);

        //输出结果
        cout << "------------------output------------------------------------" << endl;
        cout << "iteration number:" << "  " << count << endl;
        for (int i = 0; i < number; i++) {
            cout << "▲x(" << i + 1 << ") =    " << minDert[i] << endl;
        }
        cout << endl;
        cout << "the result:" << endl;
        for (int i = 0; i < number; i++) {
            cout << "x(" << i + 1 << ") =    " << nodeValue[i + 1] << endl;
        }

    }
    else if (choose == 2) {
    //------------------------同伦法求解电路方程
        double t = 0,step;
        double accurateValue;
        int count = 0, number = 0;

        cout << "Please enter the initial data number：" << endl;
        cin >> number;
        cout << "Please enter the initial data value:" << endl;
        for (int i = 0; i < number; i++) {
            cin >> nodeValue[i + 1];
            a_value[i + 1] = nodeValue[i + 1];
        }
        cout << "please input required accuracy:" << endl;
        cin >> accurateValue;

        cout << "please input step：" << endl;
        cin >> step;

        while (t <=1.0) {
            Homotopy_NR(jacMat, result, minDert, number, count, accurateValue, datum, lastnode, t);
            cout << "iteration number:" << "  " << count << endl;
            count = 0;
            t = t + step;
        }


        //输出结果
        cout << "------------------output------------------------------------" << endl;
        cout << "iteration number:" << "  " << count << endl;
        for (int i = 0; i < number; i++) {
            cout << "▲x(" << i + 1 << ") =    " << minDert[i] << endl;
        }
        cout << endl;
        cout << "the result:" << endl;
        for (int i = 0; i < number; i++) {
            cout << "x(" << i + 1 << ") =    " << nodeValue[i + 1] << endl;
        }

    }
    else if (choose == 3) {

    // -------------------------------Tran---------------------------------------------------


    outfile.open("tran.txt", ios::out);
    isTran = 1;              //--------------------
    int number = 0;

    cout << "Please enter the initial data number：" << endl;
    cin >> number;
    cout << "Please enter the initial data value:" << endl;
    for (int i = 0; i < number; i++) {
        cin >> nodeValue[i + 1];
    }

    int count = 1;
    double accurateValue;

    cout << "please input required accuracy:" << endl;
    cin >> accurateValue;
    for (double i = stepSize; i < stopTime + stepSize; i = i + stepSize) {

        for (int i = 0; i < number; i++) {
            for (int j = 0; j < number; j++) {
                jacMat[i + 1][j + 1] = 0.0;
            }
            result[i + 1] = 0.0;
        }
        stepNum++;
        nodePtr = nodeList.getNode(0);
        while (nodePtr != NULL) {
            if (nodePtr->getNameNum() != datum) {
                nodePtr->printNodalMat(datum, lastnode, result);
            }
            nodePtr = nodePtr->getNext();
        }

        compPtr = compList.getComp(0);
        while (compPtr != NULL) {
            compPtr->specialPrintMat(datum, result);
            compPtr = compPtr->getNext();
        }


        //~> go down the component list and give supernode equations for all float sources (Nodal Analysis)
        if (eqType != Modified) {
            compPtr = compList.getComp(0);
            while (compPtr != NULL) {
                compPtr->printSuperNodeMat(datum, lastnode, result);
                compPtr = compPtr->getNext();
            }
        }


        // go down the node list and give additional MNA equations
        if (eqType == Modified) {
            nodePtr = nodeList.getNode(0);
            while (nodePtr != NULL) {
                if (nodePtr->getNameNum() != datum)
                    nodePtr->printMNAMat(datum, lastnode, result);
                nodePtr = nodePtr->getNext();
            }
        }

        //求jac矩阵

        nodePtr1 = nodeList.getNode(0);
        while (nodePtr1 != NULL) {
            if (nodePtr1->getNameNum() != datum) {
                nodePtr2 = nodeList.getNode(0);
                while (nodePtr2 != NULL) {
                    if (nodePtr2->getNameNum() != datum) {
                        nodePtr1->printJacMat(datum, nodePtr2, lastnode, eqType, jacMat);
                    }
                    nodePtr2 = nodePtr2->getNext();
                }
            }
            nodePtr1 = nodePtr1->getNext();
        }

        // go down the component list and give equations for all sources
        compPtr = compList.getComp(0);
        while (compPtr != NULL) {
            nodePtr2 = nodeList.getNode(0);
            compPtr2 = compList.getComp(0);
            while (nodePtr2 != NULL) {
                if (nodePtr2->getNameNum() != datum) {
                    compPtr->specialPrintJacMat(datum, nodePtr2, lastnode, eqType, compPtr2, &specPrintJacMNA, jacMat); // ~> specPrintJacMNA is used to verify if the jacobians w.r.t. the Modified equations was already printed to print only once.
                }
                nodePtr2 = nodePtr2->getNext();
            }
            specPrintJacMNA = 0;
            compPtr = compPtr->getNext();
        }




        // print the Jacobians for the additional MNA equations
        if (eqType == Modified) {
            nodePtr1 = nodeList.getNode(0);
            while (nodePtr1 != NULL) {
                if (nodePtr1->getNameNum() != datum) {
                    nodePtr2 = nodeList.getNode(0);
                    while (nodePtr2 != NULL) {
                        if (nodePtr2->getNameNum() != datum)
                            nodePtr1->printJacMNAMat(datum, nodePtr2, lastnode, jacMat);
                        nodePtr2 = nodePtr2->getNext();
                    }
                }
                nodePtr1 = nodePtr1->getNext();
            }
        }

        //向后欧拉法
        for (int j = 1; j <= number; j++) {
            preX[j] = nodeValue[j];
        }

        NR_Iterations(jacMat, result, minDert, number, count, accurateValue, datum, lastnode);

        for (int j = 1; j <= number; j++) {
            outfile << "Time:" << i << ",x" << j << "=" << nodeValue[j] << endl;
        }
        outfile << "----------------------------" << endl;
    }

    outfile.close();

    }
    else {
    cout << "input error!!!" << endl;
    }



    return 0;
}


void NR_Iterations(double jacMat[][30], double result[], double minDert[], int number, int& count, double accurateValue, int datum
    , int lastnode, bool Homotopy, double t) {

    VectorXd F(number);
    MatrixXd Jac(number, number);
    VectorXd delta(number);
    
    get_circuit_equation(datum, lastnode);
    get_Jacobian_matrix(datum, lastnode);

    for (int i = 0; i < number; i++) {
        F(i) = result[i + 1];
    }
    for (int i = 0; i < number; i++) {
        for (int j = 0; j < number; j++) {
            Jac(i, j) = jacMat[i + 1][j + 1];
        }
    }
    

    for (int i = 0; i < number; i++) {
        F(i) = result[i + 1];
    }
    for (int i = 0; i < number; i++) {
        for (int j = 0; j < number; j++) {
            Jac(i, j) = jacMat[i + 1][j + 1];
        }
    }

    delta = Jac.fullPivLu().solve(-F);
    cout << "xxxxxxxxxxxxxxxxxxxx" << endl;
    for (int i = 0; i < number; i++) {
        minDert[i] = delta(i);
        cout <<"当前误差为：" << delta(i) << endl;
    }


    cout << "xxxxxxxxxxxxxxxxxxxx" << endl;
    for (int i = 0; i < number; i++) {
        nodeValue[i + 1] = nodeValue[i + 1] + minDert[i];
        cout << "当前值为：" << nodeValue[i + 1] << endl;
    }//更新X

    for (int i = 0; i < number; i++) {
        for (int j = 0; j < number; j++) {
            jacMat[i + 1][j + 1] = 0.0;
        }
        result[i + 1] = 0.0;
    }

    while (!isAccurate(minDert, number, accurateValue)) {

        count++;
        get_circuit_equation(datum, lastnode);
        get_Jacobian_matrix(datum, lastnode);

        for (int i = 0; i < number; i++) {
            F(i) = result[i + 1];
        }
        for (int i = 0; i < number; i++) {
            for (int j = 0; j < number; j++) {
                Jac(i, j) = jacMat[i + 1][j + 1];
            }
        }
        delta = Jac.fullPivLu().solve(-F);
        for (int i = 0; i < number; i++) {
            minDert[i] = delta(i);
        }

        for (int i = 0; i < number; i++) {
            nodeValue[i + 1] = nodeValue[i + 1] + minDert[i];
        }

        for (int i = 0; i < number; i++) {
            for (int j = 0; j < number; j++) {
                jacMat[i + 1][j + 1] = 0.0;
            }
            result[i + 1] = 0.0;
        }
    }
}


void Homotopy_NR(double jacMat[][30], double result[], double minDert[], int number, int& count, double accurateValue, int datum
    , int lastnode, double t) {

    VectorXd F(number);
    MatrixXd Jac(number, number);
    VectorXd delta(number);

    for (int i = 0; i < number; i++) {
        for (int j = 0; j < number; j++) {
            jacMat[i + 1][j + 1] = 0.0;
        }
        result[i + 1] = 0.0;
    }
    get_circuit_equation(datum, lastnode);
    get_Jacobian_matrix(datum, lastnode);
    for (int i = 1; i <= number; i++) {
        result[i] = t * result[i] + (1 - t) * 1e-3 * (nodeValue[i] - a_value[i]);
    }
    for (int i = 0; i < number; i++) {
        for (int j = 0; j < number; j++) {
            if (i == j) {
                jacMat[i + 1][j + 1] = t * jacMat[i + 1][j + 1] + (1 - t) * 1e-3;
            }
            else {
                jacMat[i + 1][j + 1] = jacMat[i + 1][j + 1] * t;
            }
        }
    }

    for (int i = 0; i < number; i++) {
        F(i) = result[i + 1];
    }
    for (int i = 0; i < number; i++) {
        for (int j = 0; j < number; j++) {
            Jac(i, j) = jacMat[i + 1][j + 1];
        }
    }

    delta = Jac.fullPivLu().solve(-F);
    //cout << "xxxxxxxxxxxxxxxxxxxx" << endl;
    for (int i = 0; i < number; i++) {
        minDert[i] = delta(i);
        //cout << "当前误差为：" << delta(i) << endl;
    }


    //cout << "xxxxxxxxxxxxxxxxxxxx" << endl;
    for (int i = 0; i < number; i++) {
        nodeValue[i + 1] = nodeValue[i + 1] + minDert[i];
        //cout << "当前值为：" << nodeValue[i + 1] << endl;
    }//更新X
    count++;

    while (!isAccurate(minDert, number, accurateValue)) {

        count++;
        for (int i = 0; i < number; i++) {
            for (int j = 0; j < number; j++) {
                jacMat[i + 1][j + 1] = 0.0;
            }
            result[i + 1] = 0.0;
        }
        get_circuit_equation(datum, lastnode);
        get_Jacobian_matrix(datum, lastnode);
        for (int i = 1; i <= number; i++) {
            result[i] = t * result[i] + (1 - t) * 1e-3 * (nodeValue[i] - a_value[i]);
        }
        for (int i = 0; i < number; i++) {
            for (int j = 0; j < number; j++) {
                if (i == j) {
                    jacMat[i + 1][j + 1] = t * jacMat[i + 1][j + 1] + (1 - t) * 1e-3;
                }
                else {
                    jacMat[i + 1][j + 1] = t * jacMat[i + 1][j + 1];
                }
            }
        }


        for (int i = 0; i < number; i++) {
            F(i) = result[i + 1];
        }
        for (int i = 0; i < number; i++) {
            for (int j = 0; j < number; j++) {
                Jac(i, j) = jacMat[i + 1][j + 1];
            }
        }
        delta = Jac.fullPivLu().solve(-F);
        for (int i = 0; i < number; i++) {
            minDert[i] = delta(i);
        }

        for (int i = 0; i < number; i++) {
            nodeValue[i + 1] = nodeValue[i + 1] + minDert[i];
        }

    }
}


bool isAccurate(double minDert[], int num, double acc) {
    bool re = true;
    for (int i = 0; i < num; i++) {
        if (minDert[i] > acc || -minDert[i] > acc) {
            re = false;
        }
    }
    return re;

}
//求电路方程
void get_circuit_equation(int datum, int lastnode) {
    Component* compPtr, * compPtr2;
    Node* nodePtr, * nodePtr1, * nodePtr2;
    int specPrintJacMNA = 0;
    EquaType eqType = Modified;

    nodePtr = nodeList.getNode(0);
    while (nodePtr != NULL) {
        if (nodePtr->getNameNum() != datum) {
            nodePtr->printNodalMat(datum, lastnode, result);
        }
        nodePtr = nodePtr->getNext();
    }

    compPtr = compList.getComp(0);
    while (compPtr != NULL) {
        compPtr->specialPrintMat(datum, result);
        compPtr = compPtr->getNext();
    }


    //~> go down the component list and give supernode equations for all float sources (Nodal Analysis)
    if (eqType != Modified) {
        compPtr = compList.getComp(0);
        while (compPtr != NULL) {
            compPtr->printSuperNodeMat(datum, lastnode, result);
            compPtr = compPtr->getNext();
        }
    }


    // go down the node list and give additional MNA equations
    if (eqType == Modified) {
        nodePtr = nodeList.getNode(0);
        while (nodePtr != NULL) {
            if (nodePtr->getNameNum() != datum)
                nodePtr->printMNAMat(datum, lastnode, result);
            nodePtr = nodePtr->getNext();
        }
    }
}

//求雅各比矩阵
void get_Jacobian_matrix(int datum, int lastnode) {
    Component* compPtr, * compPtr2;
    Node* nodePtr, * nodePtr1, * nodePtr2;
    int specPrintJacMNA = 0;
    EquaType eqType = Modified;

    nodePtr1 = nodeList.getNode(0);
    while (nodePtr1 != NULL) {
        if (nodePtr1->getNameNum() != datum) {
            nodePtr2 = nodeList.getNode(0);
            while (nodePtr2 != NULL) {
                if (nodePtr2->getNameNum() != datum) {
                    nodePtr1->printJacMat(datum, nodePtr2, lastnode, eqType, jacMat);
                }
                nodePtr2 = nodePtr2->getNext();
            }
        }
        nodePtr1 = nodePtr1->getNext();
    }

    // go down the component list and give equations for all sources
    compPtr = compList.getComp(0);
    while (compPtr != NULL) {
        nodePtr2 = nodeList.getNode(0);
        compPtr2 = compList.getComp(0);
        while (nodePtr2 != NULL) {
            if (nodePtr2->getNameNum() != datum) {
                compPtr->specialPrintJacMat(datum, nodePtr2, lastnode, eqType, compPtr2, &specPrintJacMNA, jacMat); // ~> specPrintJacMNA is used to verify if the jacobians w.r.t. the Modified equations was already printed to print only once.
            }
            nodePtr2 = nodePtr2->getNext();
        }
        specPrintJacMNA = 0;
        compPtr = compPtr->getNext();
    }




    // print the Jacobians for the additional MNA equations
    if (eqType == Modified) {
        nodePtr1 = nodeList.getNode(0);
        while (nodePtr1 != NULL) {
            if (nodePtr1->getNameNum() != datum) {
                nodePtr2 = nodeList.getNode(0);
                while (nodePtr2 != NULL) {
                    if (nodePtr2->getNameNum() != datum)
                        nodePtr1->printJacMNAMat(datum, nodePtr2, lastnode, jacMat);
                    nodePtr2 = nodePtr2->getNext();
                }
            }
            nodePtr1 = nodePtr1->getNext();
        }
    }
}

double stripString(char* stringIn) {
    char buf[BufLength], buf2[BufLength];
    int a, b;
    strcpy(buf, stringIn);
    for (a = 0; buf[a] != '='; a++) {};
    a++;
    for (b = 0; buf[a] != '\0'; b++, a++)
        buf2[b] = buf[a];
    buf2[b] = '\0';
    return atof(buf2);    //atof()函数 将字符转换为浮点数型
};

//Print the linked list of components to check
void printComponents(Component* compPtr) {
    char compTypeName[6];
    cout << endl << "Components: " << endl << endl;
    while (compPtr != NULL) {
        strcpy(compTypeName, strComponentType(compPtr));
        cout << "->" << compTypeName << compPtr->getcompNum();
        compPtr = compPtr->getNext();
    }
    cout << endl;
    return;
}


void printNodes(Node* nodePtr, int compFlag) {

    Connections* conPtr;
    cout << endl << "Nodes: " << endl << endl;
    while (nodePtr != NULL) {
        if (compFlag == 0) { //It is printed just the names of the nodes
            cout << "-> " << nodePtr->getNameNum();
        }
        else if (compFlag == 1) { //It is printed the nodes and the connections
            cout << "-> " << nodePtr->getNameNum() << " {";
            conPtr = nodePtr->getConList();
            while (conPtr->next != NULL) {
                cout << strComponentType(conPtr->comp) << conPtr->comp->getcompNum() << ", ";
                conPtr = conPtr->next;
            }
            cout << strComponentType(conPtr->comp) << conPtr->comp->getcompNum() << '}' << endl;
        }
        else {
            cout << "Invalid value for compFlag. (0) to print just nodes, (1) to print nodes and connections!";
            exit(1);

        }

        nodePtr = nodePtr->getNext();
    }


    return;
}


char* strComponentType(Component* compPtr) {

    char* compTypeName = new char[6];
    switch (compPtr->getType()) {

    case VSource: strcpy(compTypeName, "V"); break;
    case Resistor: strcpy(compTypeName, "R"); break;
    case BJT: strcpy(compTypeName, "T"); break;
    case MOSFET: strcpy(compTypeName, "M"); break;
    case ISource: strcpy(compTypeName, "I"); break;
    case Inductor: strcpy(compTypeName, "ind"); break;
    case Diode: strcpy(compTypeName, "Diode"); break;
    case Capacitor: strcpy(compTypeName, "C"); break;
    }

    return compTypeName;
}



char* ComponentTypeName(Component* compPtr) {

    char* compTypeName = new char[6];
    switch (compPtr->getType()) {

    case VSource: strcpy(compTypeName, "VSource"); break;
    case Resistor: strcpy(compTypeName, "Resistor"); break;
    case BJT: strcpy(compTypeName, "BJT"); break;
    case MOSFET: strcpy(compTypeName, "MOSFET"); break;
    case ISource: strcpy(compTypeName, "ISource"); break;
    case Inductor: strcpy(compTypeName, "Inductor"); break;
    case Diode: strcpy(compTypeName, "Diode"); break;
    case Capacitor: strcpy(compTypeName, "Capacitor"); break;
    }

    return compTypeName;
}


int portNum(Component* comPtr, Node* nodePtr) {
    if (comPtr->getNodeNum(0) == nodePtr->getNum()) {
        return 0;
    }
    else {
        return 1;
    }
}

void LU(double A[][30], double x[], double b[], int n) {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrixA(n, n);
    Eigen::Matrix<double, Eigen::Dynamic, 1> vectorb(n);
    Eigen::PartialPivLU<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> lu_decomp(n);

    // 将数组 A 的值赋给矩阵 matrixA
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrixA(i, j) = A[i][j];
        }
    }

    // 将数组 b 的值赋给向量 vectorb
    for (int i = 0; i < n; i++) {
        vectorb(i) = b[i];
    }

    lu_decomp.compute(matrixA);  // 对矩阵进行 LU 分解
    Eigen::Matrix<double, Eigen::Dynamic, 1> x_dert = lu_decomp.solve(vectorb);  // 求解方程 Ax=b

    // 将结果复制到数组 x
    for (int i = 0; i < n; i++) {
        x[i] = x_dert(i);
    }
}

void convertArray(double jacMat[][30], double A[][30], double result[], double y[], int number) {
    for (int i = 0; i < number; i++) {
        for (int j = 0; j < number; j++) {
            A[i][j] = jacMat[i + 1][j + 1];
        }

        y[i] = -result[i + 1];
    }
}



















