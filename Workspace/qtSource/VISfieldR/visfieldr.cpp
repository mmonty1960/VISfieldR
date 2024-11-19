/*Main author: Marco Montecchi
             ENEA (Italy)
             email: marco.montecchi@enea.it
Porting to Windows and advanced oscillators by
             Alberto Mittiga
             ENEA (Italy)
             email: alberto.mittiga@enea.it


VISfieldR is the rotational variant of the VISfield method,
designed for the evaluation of the intercept-factor of parabolic-trough
modules from a series of images captured from the same point of view.
This viewpoint is situated on the plane that intersects the module's
rotation axes and contains an horizontal line.
In this rotational version, the image sequence is obtained as the module
rotates in such a way that the image of the receiver tube scan the
parabolic surface from one rim to the other.



   Copyright (C) 2025  Marco Montecchi

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#include "visfieldr.h"
#include <QtGui>
#include <QFile>
#include <QFileDialog>
#include <QMessageBox>
#include <QKeyEvent>
#include <fstream>
#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/calib3d.hpp>
#include <cminpack.h>
#include <QFontDialog>
#include <QFontInfo>

using namespace std;
using namespace cv;

//global variables *******************************
QString pathbase;//="/run/media/marco/d2quadra/VISfieldR/";//"/home/marco/Workspace/VISfieldR/";
QString fileInfo="";
QString fileResult=pathbase+"results.txt";
QString pathfile;
QString PanelType;
QString msgPoint[9];
QString bNf,eImg,msg,imgNameVertex="";
int Nraw,N_Phang;
int mpxX1,mpxY1,mpxX2,mpxY2,drag,soglia;
int jROImin,iROImin,Njrs,Nirs;
int pippo,iPanel,ioBrek;
int Nt=8;//Ntargets= 2joint-pivot + 4corners +2receiver
int Np=0;//Npoints
int Nmin,Nmax;//N frame Min Max
int iFull0Half1=0;
double px[99][4]; //coordinate in pixel dei target nell'immagine e del punto calcolato
double P[99][3];//coordinate x,y,z targets
double focal,panel_width,xmin,xmax,length,lengthUP,lengthDW,gapUP,gapDW,radius,df=0.;
double xrMin,xrMax,Magne,xsMin,xsMax,xP,zP,dview,xStep,xCa,yCa,zCa;
double Phold[12][2];
double corrPhold[20][12];//[iz][ix]=correction to Phold at the croos of iz panel-raw and ix column of hanging point
double matAlpha[20][4];//[iz][ix]=Alpha matrix at the croos of iz-raw of z and ix panel-column of x
double chi2,chi2fin;
double pig=3.141592654;
double Vservice[3];//vettore di servizio
double thetaImg;
double IntFact[600][2];
double mjt,mjb,mil,mir,cjt,cjb,cil,cir,dpximg,iLeft,jLeft,jRight,iRight,Lpx;
// reference frame in the parabolic-trough module
//    x ^
//      |
//      |
//      O----->y


// Nikon 1 parameters **************************************************************************
int Width= 5232;
int Height=3488;
int NpxX=5232;//N. px lato X                   (0,0) ---> x' (j)
int NpxY=3488;//N. px lato Y                     |
double ccdY=8.8   ;//lunghezza lato Y ccd        |
double ccdX=13.2  ;//lunghezza lato X ccd        V y' (i)
//calib dic 2017 con target F81 ********************************
//double fcam=10.290;//10.293; //camera focal (mm)
//double xp_iWitness=0.012;//0.010;
//double yp_iWitness=-0.189;//-0.202;
//double k1_iWitness=1.528E-03;//1.489E-03;
//double k2_iWitness=-1.159E-05;//-8.620E-06;
//double k3_iWitness=8.615E-08;//1.456E-08;
//double p1_iWitness=2.287E-04;//2.041E-04;
//double p2_iWitness=3.442E-05;//6.411E-05;
//calib 25/gen/2018 con target su modulo PCS********************
double fcam=10.2849;// camera focal (mm)
double xp_iWitness=0.0157;
double yp_iWitness=-0.1778;
double k1_iWitness=1.5428E-03;
double k2_iWitness=-1.0657E-05;
double k3_iWitness=3.2066E-08;
double p1_iWitness=1.8961E-04;
double p2_iWitness=1.1985E-04;
//***************************************************************
double pxdimX=ccdX/double(NpxX);
double pxdimY=ccdY/double(NpxY);
double fx=fcam/pxdimX;
double fy=fcam/pxdimY;
double cx=NpxX/2.+xp_iWitness/pxdimX;
double cy=NpxY/2.-yp_iWitness/pxdimY;

Mat cameraMatrix =  (Mat_<double>(3,3) << fcam/pxdimX,          0 ,     cx,
                                             0,       fcam/pxdimY ,     cy,
                                             0,                 0 ,      1);

Mat distCoeffs = (Mat_<double>(5,1)   <<    -k1_iWitness*pow(fcam,2.),-k2_iWitness*pow(fcam,4.),
                  +p1_iWitness*pow(fcam,1.),-p2_iWitness*pow(fcam,1.),-k3_iWitness*pow(fcam,6.));
//***************************************************************************************************
Point2i Pt1,Pt2;
Mat img;// image
Mat Cimg;//cropped image
Mat imgDisplayed;//image displayed
Mat smap, smapSlope;
Mat img8b(Height,Width,CV_8UC1);
Mat img8RGB(Height,Width,CV_8UC3);
int NiMappa=500;
int NjMappa=1000;
int arr[3] = {NiMappa,NjMappa,6};
    Mat mappa(3, arr, CV_64FC1);
    //[nf][j][0]=counter
    //[nf][j][1]=gray level
    //[nf][j][2]=dev_arctangent
    //[nf][j][3]=Int Fat
    //[nf][j][4]=cos(theta)
    //[nf][j][5]=sensitivity
int cropInfo[5];
double pxX,pxY;//coordinate pixel i,j
int irsStep=1;
int irsStart=0;
int iProcess=0;
Vec3b vRGB;

// camera position & attitude
double Pc[3];//coordinate camera drone
double EA[3];//yaw pitch roll camera
double xSto[6];//position & attitude at reference frame
int Nimg,NimgRef,iSetRef=0;
double valT,valT0;

//invoked functions ********************************
int fcn(void *p, int m, int n, const double *x, double *fvec, int iflag);
void xyz2pxpy(double X, double Y, double Z);
void Mod2Dro(double yaw,double pitch, double roll,double xm, double ym, double zm);
double zero(int ic, double tetaMIN, double tetaMAX, double xr);
void imgDimension(double tg2p);
int segno(double val);
void setRGB(double x);

struct pointToFit{
 int Nt;
}pTF[1];



void on_mouse(int ev, int x, int y, int, void* obj){
  VISfieldR* app = static_cast<VISfieldR*>(obj);
  if (app)
    app->on_mouse_internal(ev, x, y);
}


VISfieldR::VISfieldR(QWidget *parent){
    setupUi(this);
    this-> setCursor(Qt::ArrowCursor);

    // signals/slots mechanism in action
    connect( FileInfo,     SIGNAL(textChanged(QString)),  this, SLOT( setPanel()));
    connect( selectFile,   SIGNAL( clicked() ),           this, SLOT( SetFile()));
    connect( selectFile_2, SIGNAL( clicked() ),           this, SLOT( SetMeasFile()));
    connect( QDSBlengthUP, SIGNAL( valueChanged(double)), this, SLOT( setGap()));
    connect( QDSBlengthDW, SIGNAL( valueChanged(double)), this, SLOT( setGap()));
    connect(DSBslopeThres, SIGNAL( valueChanged(double)), this, SLOT( plotMap()));
    connect( sB_Nimg ,     SIGNAL( valueChanged(int)),    this, SLOT( viewFrameX()));
    connect( sB_gLup,      SIGNAL( valueChanged(int)),    this, SLOT( viewFrameX()));
    connect( sB_gLdown,    SIGNAL( valueChanged(int)),    this, SLOT( viewFrameX()));
    connect( checkBoxGray, SIGNAL( stateChanged(int)),    this, SLOT( viewFrameX()));
    connect( pB_setCamera, SIGNAL( clicked() ),           this, SLOT( setCamera()));
    connect( pB_crop,      SIGNAL( clicked() ),           this, SLOT( crop()));
    connect( pB_crop_IF,   SIGNAL( clicked() ),           this, SLOT( cropIF()));
    connect( pB_crop_devSlope,SIGNAL( clicked() ),        this, SLOT( cropDevSlope()));
    connect( pB_setJPp,    SIGNAL( clicked() ),           this, SLOT( setJPp()));
    connect( pB_saveSet,   SIGNAL( clicked() ),           this, SLOT( saveMeasFile()));
    connect( pB_process,   SIGNAL( clicked() ),           this, SLOT( process()));
    connect( sB_Nmin,      SIGNAL( valueChanged(int)),    this, SLOT( updateN()));
    connect( sB_Nmax,      SIGNAL( valueChanged(int)),    this, SLOT( updateN()));
    connect( sB_NimgRef,   SIGNAL( valueChanged(int)),    this, SLOT( updateN()));
    connect( pB_viewNmin,  SIGNAL( clicked() ),           this, SLOT( viewNmin()));
    connect( pB_viewNmax,  SIGNAL( clicked() ),           this, SLOT( viewNmax()));
    connect( pB_viewImgRef,SIGNAL( clicked() ),           this, SLOT( viewImgRef()));
    connect(pushButton_setFont, SIGNAL( clicked() ), this, SLOT(setFontDia()));
    connect( doubleSpinBox_offset, SIGNAL( valueChanged(double)), this, SLOT( viewFrameX()));
    connect(comboBox_method, SIGNAL(currentIndexChanged(int)), this, SLOT(enableDisable()));

    //create OpenCV windows
    namedWindow("Image",WINDOW_NORMAL);
    namedWindow("ROI",WINDOW_NORMAL);
    namedWindow("Module map",WINDOW_NORMAL);
    namedWindow("Module Slope-Deviation map", WINDOW_NORMAL);
    imshow("Image",img8RGB);
    imshow("ROI",img8RGB);
    setWin("Image");//abilita il cursore sulla finestra
    setWin("ROI");
    setWin("Module map");
    setWin("Module Slope-Deviation map");
    resizeWindow("Image",int(Width/5),int(Height/5));
    resizeWindow("ROI",int(Width/5),int(Height/5));
    waitKey(10);

    msgPoint[0]="Select the ROI containing the joint-pivot on the RIGHT";
    msgPoint[1]="Select the ROI containing the joint-pivot on the LEFT";
    msgPoint[2]="select the ROI containing the corner DOWN RIGHT";
    msgPoint[3]="select the ROI containing the corner TOP RIGHT";
    msgPoint[4]="select the ROI containing the corner TOP LEFT";
    msgPoint[5]="select the ROI containing the corner DOWN LEFT";
    msgPoint[6]="select the ROI containing the receiver RIGHT";
    msgPoint[7]="select the ROI containing the receiver LEFT";
    msgPoint[8]="Done!";

    ioBrek=0;
    printf("cx=%f cy=%f\n",cx,cy);
    setFontDia();
    enableDisable();
    pB_crop_IF -> setEnabled(false);
    pB_crop_devSlope -> setEnabled(false);

    // parameter initialization
#ifdef __unix__
#define IS_POSIX 1
#else
#define IS_POSIX 0
#endif
    QDir dir2;  //current directory
    dir2.cdUp();//cd ..
    dir2.cdUp();//cd ..
    if (IS_POSIX == 1) {
        //Linux path initialization
        //nothing to do
    }
    else {
        //windows path inizialization
        dir2.cdUp();//cd ..
    }
    pathbase=dir2.absolutePath()+"/VISfieldR/";
    printf("dir= %s\n",pathbase.toStdString().c_str());
}

VISfieldR::~VISfieldR()
{
    delete ui;
}

void VISfieldR::enableDisable(){
    int iMeth=comboBox_method->currentIndex();
    if(iMeth==0){
        doubleSpinBox_offset->setEnabled(true);
        dSB_thetaStep->setEnabled(true);
        dSB_mt->setEnabled(false);
        dSB_q->setEnabled(false);
        lineEdit_t1->setEnabled(false);
    }
    else{
        doubleSpinBox_offset->setEnabled(false);
        dSB_thetaStep->setEnabled(false);
        dSB_mt->setEnabled(true);
        dSB_q->setEnabled(true);
        lineEdit_t1->setEnabled(true);
    }
}

void VISfieldR::setFontDia(){
    bool ok=true;
    //QFont myFont=QFontDialog::getFont(&ok,this);//,QFont("Noto Sans",8)
    QFont myFont("Comic Sans MS", 14);
    if(ok) {
        // the user clicked OK and font is set to the font the user selected
        const QWidgetList allWidgets = QApplication::allWidgets();
        for (QWidget *widget : allWidgets){
            widget->setFont(myFont);
            widget->update();
        }
        QCoreApplication::processEvents();
    } else {
        // the user canceled the dialog; font is set to the initial
        // value, in this case Helvetica [Cronyx], 10
    }
}


void VISfieldR::setWin(const std::string& _winname)
{
  cv::namedWindow(_winname);
  this-> winname = _winname;
  cv::setMouseCallback(winname, on_mouse, this);
}



void VISfieldR::on_mouse_internal(int event, int x, int y){
    Point point;
    //uint16_t value=0;
    int value;
    if(x>=0 && x<Width && y>=0 && y<Height){
        lineEdit_j -> setText(QString::number(x)+" -> "+QString::number(x-cx));
        lineEdit_i -> setText(QString::number(y)+" -> "+QString::number(y-cy));
        //value=img.at<uint16_t>(y,x);
        value=img8b.at<uchar>(y,x);
        lineEdit_GL -> setText(QString::number(value));
    }
    if (event == EVENT_LBUTTONDOWN){
        point = Point(x, y);
        if(x<0) x=0;
        if(x>Width) x=Width;
        if(y<0) y=0;
        if(y>Height) y=Height;
        mpxX1=x;//j
        mpxY1=y;//i
        drag=1;
        //printf("left button DW j=%d i=%d\n",mpxX1,mpxY1);
    }
    if (event == EVENT_LBUTTONUP){
        point = Point(x, y);
        if(x<0) x=0;
        if(x>Width) x=Width;
        if(y<0) y=0;
        if(y>Height) y=Height;
        mpxX2=x;
        mpxY2=y;
        if(mpxX2 > mpxX1 && mpxY2 > mpxY1)
         drag=2;
        else
         drag=1;
        //printf("left button UP j=%d i=%d\n",mpxX2,mpxY2);
    }
    if (event == EVENT_RBUTTONDOWN) {
      drag=-1;
      point = Point(x, y);
      if(x<0) x=0;
      if(x>Width) x=Width;
      if(y<0) y=0;
      if(y>Height) y=Height;
      mpxX1=x;
      mpxY1=y;
      //value=img.at<uint16_t>(mpxY1,mpxX1);
      value=img8b.at<uchar>(mpxY1,mpxX1);
      //printf("right button DW j=%d i=%d grayLevel=%d\n",mpxX1,mpxY1,value);
    }
    if (event == EVENT_MBUTTONDOWN) {
        drag=10;
        point = Point(x, y);
        mpxX1=x+jROImin;
        mpxY1=y+iROImin;
        sB_j->setValue(mpxX1);
        sB_i->setValue(mpxY1);
        //printf("(i,j) = (%d , %d)\n",mpxY1,mpxX1);
        fflush(stdout);
        if(mpxY1<Height && mpxX1<Width){
            //value=img.at<uint16_t>(mpxY1,mpxX1);
            value=img8b.at<uchar>(mpxY1,mpxX1);
            lineEdit_j -> setText(QString::number(mpxX1));
            lineEdit_i -> setText(QString::number(mpxY1));
            lineEdit_GL-> setText(QString::number(value));
            //printf("medium button DW j=%d i=%d grayLevel=%d\n",mpxX1,mpxY1,value);
        }else
            drag=0;
    }
}


void VISfieldR::SetFile(){
    QString path2 = QFileDialog::getOpenFileName(this,tr("Choose the InfoFile"),
        pathbase + "InfoFile");
     if(!path2.isEmpty()){
      QString name = QFileInfo(path2).fileName();
      FileInfo -> setText( name );
     }
}


//settaggio pannello
void VISfieldR::setPanel(){
    if(ioBrek==1)
        return;
    else
        ioBrek=1;
    QString stringa,fname,fileInfo,comment;
    double dval,pagx,pagy;
    printf("->setPanel()...\n");
    fname = FileInfo -> text();
    fileInfo=pathbase + "InfoFile/"+fname;
    printf("setPanel by file=%s\n",fileInfo.toStdString().c_str());
    QFile fileP(fileInfo);
    QTextStream fileA(&fileP);
    if(!fileP.open(QIODevice::ReadOnly | QIODevice::Text)){
        msg="Can not load "+fileInfo;
        textEdit -> append (msg);
        return;
    }
      fileA >> comment >>  stringa;
      label_manufacturer ->setText(stringa);
      fileA >> comment >>  stringa;
      labe_Ptype ->setText(stringa);
      fileA >> comment >>  focal;
      QDSBfocale -> setValue(focal);
      fileA >> comment >>  panel_width;
      QDSBwidth -> setValue(panel_width);
      fileA >> comment >>  dval;
      QDSBinnerXmin -> setValue(dval);
      xmin=dval;
      fileA >> comment >>  dval;
      QDSBinnerXmax -> setValue(dval);
      if(stringa.toStdString()=="P/2"){
       iPanel=0;
       PanelType="P/2";
       xmax=dval;
       QDSBouterXmin -> setEnabled ( false );
       QDSBouterXmax -> setEnabled ( false );
       QDSBmiddleXmin -> setEnabled ( false );
       QDSBmiddleXmax -> setEnabled ( false );
       QDSBouterXmin -> setValue(0.0);
       QDSBouterXmax -> setValue(0.0);
       QDSBmiddleXmin -> setValue(0.0);
       QDSBmiddleXmax -> setValue(0.0);
      }
      else if(stringa.toStdString()=="P/4"){
       iPanel=1;
       PanelType="P/4";
       QDSBouterXmin -> setEnabled ( true );
       QDSBouterXmax -> setEnabled ( true );
       QDSBmiddleXmin -> setEnabled ( false );
       QDSBmiddleXmax -> setEnabled ( false );
       QDSBmiddleXmin -> setValue(0.0);
       QDSBmiddleXmax -> setValue(0.0);
       fileA >> comment >>  dval;
       QDSBouterXmin -> setValue(dval);
       fileA >> comment >>  dval;
       QDSBouterXmax -> setValue(dval);
       xmax=dval;
      }
      else if(stringa.toStdString()=="P/6"){
       iPanel=2;
       PanelType="P/6";
       QDSBouterXmin -> setEnabled ( true );
       QDSBouterXmax -> setEnabled ( true );
       QDSBmiddleXmin -> setEnabled ( true );
       QDSBmiddleXmax -> setEnabled ( true );
       fileA >> comment >>  dval;
       QDSBouterXmin -> setValue(dval);
       fileA >> comment >>  dval;
       QDSBouterXmax -> setValue(dval);
       xmax=dval;
       fileA >> comment >>  dval;
       QDSBmiddleXmin -> setValue(dval);
       fileA >> comment >>  dval;
       QDSBmiddleXmax -> setValue(dval);
      }
      fileA >> comment;
      N_Phang=-1;//N Phang raw for panel; a couple of hanging points is assumed to be at each raw
      while (!fileA.atEnd()) {//each loaded hanging points is supposed to be paired with another at same x and -z
          N_Phang++;
          fileA >> pagx >> pagy;
          Phold[N_Phang][0]=pagx;//x
          Phold[N_Phang][1]=pagy;//z
      }
      fileP.close();
      int idata=N_Phang;
      Nraw=0;//N raw of panels
      if(iPanel == 0){
          Nraw=1;
      }
      else if(iPanel == 1){
          //N_Phang=N_Phang/2;
          Nraw=2;
      }
      else if(iPanel == 2){
          //N_Phang=N_Phang/3;
          Nraw=3;
      }
      if(iFull0Half1==0){
          Nraw=Nraw*2;
          for(int i=idata; i<2*idata; i++){
              Phold[i][0]=Phold[i-idata][0];
              Phold[i][1]=Phold[i-idata][1];
          }
          for(int i=0; i<idata; i++){
              Phold[i][0]=-Phold[2*idata-1-i][0];
              Phold[i][1]=Phold[2*idata-1-i][1];
          }
      }
      printf("Nraw= %d N_Phang= %d\n",Nraw,N_Phang);
      for(int i=0; i<N_Phang*Nraw; i++){
          printf("Phold[%d][0]=%f\tPhold[%d][1]=%f\n",i,Phold[i][0],i,Phold[i][1]);
      }
      ioBrek=0;
      setGap();
}


void VISfieldR::setGap (){
    if(ioBrek==1)
        return;
    else
        ioBrek=1;
    printf("->setGap()...\n");
    lengthUP = QDSBlengthUP -> value();
    lengthDW = QDSBlengthDW -> value();
    length=(lengthUP+lengthDW)/2.;
    panel_width = QDSBwidth -> value();
    int Npanel=int(length/(panel_width)+0.5);//N panel along one raw
    if(Npanel>1){
        gapUP=(lengthUP-double(Npanel)*panel_width)/double(Npanel-1);
        gapDW=(lengthDW-double(Npanel)*panel_width)/double(Npanel-1);
    }
    else{
        gapUP=0.;
        gapDW=0.;
    }
    printf("lenghtUP=%f lengthDW=%f panel_width=%f\nNpanel4row=%d gapUP=%f gapDW=%f\n",
           lengthUP,lengthDW,panel_width,Npanel,gapUP,gapDW);
    if(gapUP < 0.0 || gapDW <0.0){
        QString msg="Error!\ngap-value has to be POSITIVE";
        QMessageBox msgBox;
        msgBox.setText(msg);
        msgBox.exec();
    }
    QDSBgapUP -> setValue(gapUP);//gap between panel
    QDSBgapDW -> setValue(gapDW);//gap between panel
    xStep=(2.-double(iFull0Half1))*xmax/double(Nmax-Nmin);
    Njrs=int((lengthUP+lengthDW)/2./xStep);
    ioBrek=0;
}


//set the measurement file to be reprocessed
void VISfieldR::SetMeasFile(){
    printf("->SetMeasFile()");
    QString comment,strg,line,pezzo;
    QStringList List;
    int value,ix0here;
    double dvalue,Pjp[3];
    fileInfo=QFileDialog::getOpenFileName(this,tr("Choose the infoMeasurement.txt file"),pathbase,tr("infoMeasurement (*.txt)"));
    QFile file(fileInfo);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)){
        msg="Can not load "+fileInfo;
        textEdit -> append (msg);
        return;
    }
    pathfile=fileInfo.section('/',0,-2);
    if(!pathfile.isNull()){
        if(!pathfile.endsWith("/"))
            pathfile=pathfile+"/";
        lineEdit_img -> setText(pathfile.section('/', -3));
        ioBrek=1;
        iSetRef=0;//camera has to be set
        cropInfo[0]=0;//Crop region has to be set
        Np=0;//all Points must be set
        for(int i=0;i<300;i++)
            IntFact[i][0]=0.;
        QTextStream stream ( &file );
        stream>>comment>>strg;
        plant_code -> setText(strg);
        stream>>comment>>strg;
        module_code -> setText(strg);
        stream>>comment>>strg;
        FileInfo -> setText(strg);
        stream>>comment>>dvalue;
        QDSBlengthUP -> setValue(dvalue);
        stream>>comment>>dvalue;
        QDSBlengthDW -> setValue(dvalue);
        stream>>comment>>dvalue;
        dSB_Zc -> setValue(dvalue);
        stream>>comment>>bNf;
        stream>>comment>>eImg;
        stream>>comment>>value;
        sB_gLup->setValue(value);
        stream>>comment>>value;
        sB_gLdown->setValue(value);
        stream>>comment>>value;
        sB_Nimg -> setMinimum(value);
        sB_Nmin -> setValue(value);
        Nmin=value;
        stream>>comment>>value;
        sB_Nimg -> setMaximum(value);
        sB_Nmax -> setValue(value);
        Nmax=value;
        stream>>comment>>dvalue;
        dSB_mt -> setValue(dvalue);
        if(abs(dvalue)>1.e-6)
            comboBox_method->setCurrentIndex(1);
        stream>>comment>>dvalue;
        dSB_q -> setValue(dvalue);
        stream>>comment>>dvalue;
        dSB_vjp -> setValue(dvalue);
        line=stream.readLine();
        line=stream.readLine();
        line=line.simplified();
        List =line.split(" ");
        int nV=List.count();
        //printf("nV=%d\n",nV);
        if(nV==4){
            for(int k=1;k<nV;k++){
                pezzo=List.at(k).toLocal8Bit().constData();
                Pjp[k-1]=pezzo.toDouble();
                printf("Pjp[%d]=%f\n",k-1,Pjp[k-1]);
            }
            doubleSpinBox_xR->setValue(Pjp[0]);
            doubleSpinBox_yR->setValue(Pjp[1]);
            doubleSpinBox_zR->setValue(Pjp[2]);
            line=stream.readLine();
            line=line.simplified();
            List =line.split(" ");
            for(int k=1;k<nV;k++){
                pezzo=List.at(k).toLocal8Bit().constData();
                Pjp[k-1]=pezzo.toDouble();
                printf("Pjp[%d]=%f\n",k-1,Pjp[k-1]);
            }
            doubleSpinBox_xL->setValue(Pjp[0]);
            doubleSpinBox_yL->setValue(Pjp[1]);
            doubleSpinBox_zL->setValue(Pjp[2]);
            fflush(stdout);
        }
        else{
            QMessageBox msgBox;
            msgBox.setText("Old version: please manually set x,y,z of  the joint-pivot markers!");
            msgBox.exec();

        }
        if(nV==4)
            stream>>comment>>strg;
        else
            strg=List.at(1).toLocal8Bit();
        //cout<<"strg= "<<strg.toStdString()<<"\n";
        fflush(stdout);
        if(strg.contains("half", Qt::CaseInsensitive))
            iFull0Half1=1;
        else
            iFull0Half1=0;
        stream>>comment>>strg;
        if(strg.contains("no", Qt::CaseInsensitive))
            checkBox_fullLength->setChecked(false);
        else
            checkBox_fullLength->setChecked(true);
        stream>>comment>>strg;
        if(strg.contains("no", Qt::CaseInsensitive))
            checkBox_useOldPoints->setChecked(false);
        else
            checkBox_useOldPoints->setChecked(true);
        stream>>comment>>strg;
        if(strg.contains("no", Qt::CaseInsensitive)){
            checkBox_x0here->setChecked(false);
            ix0here=0;
        }else{
            checkBox_x0here->setChecked(true);
            ix0here=1;
        }
        stream>>comment>>strg;
        if(!comment.isEmpty()){
            if(ix0here==0){
                QString last=strg.section('/',-1);
                value=last.toInt();
                printf("last= %s\n",last.toStdString().c_str());
                if (value>=0){
                    NimgRef=value;
                    Np=8;
                    sB_NimgRef ->setValue(NimgRef);
                    checkBox_useOldPoints->setCheckState(Qt::Checked);
                    sB_Nimg->setValue(NimgRef);
                }else{
                    Np=0;
                }
            }else{
                imgNameVertex=strg;
                if(!imgNameVertex.isEmpty()){
                    Np=6;
                    QString last=strg.section('_',-1);
                    strg=last.section('.',0,0);
                    int iHead=bNf.size();
                    strg=strg.remove(0,iHead-1);
                    value=strg.toInt();
                    //printf("strg=%s value= %d\n",strg.toStdString().c_str(),value);
                    lineEdit_img ->setText(imgNameVertex.section('/', -4));
                    checkBox_useOldPoints->setCheckState(Qt::Checked);
                    sB_NimgRef ->setValue(value);
                    NimgRef=value;
                    if(NimgRef>=Nmin && NimgRef<=Nmax)
                        Nimg=NimgRef;
                    else
                        Nimg=int((Nmax+Nmin)/2);
                    sB_Nimg->setValue(Nimg);
                }else{
                    Np=0;
                }
            }
        }
        if(Np>0){
            for(int i=0; i<Np; i++){
                stream >> comment >> px[i][0]>>px[i][1];
                printf("px[%d][0]= %f\tpx[%d][1]= %f\n",i,px[i][0],i,px[i][1]);
            }
            stream>>comment>>dvalue;
            chi2=dvalue;
            lineEdit_chi2->setText(QString::number(chi2,'f',1));
        }
        else{
            checkBox_useOldPoints -> setCheckState(Qt::Unchecked);
            sB_NimgRef ->setValue(-1);
        }
        file.close();
        printf("The stored conspicuous points are Np=%d\n",Np);
        ioBrek=0;
        setPanel();

        //show maps if already computed
        QString mapImg=pathfile+"module_map.jpeg";
        smap=imread(mapImg.toStdString(),IMREAD_UNCHANGED);
        if(!smap.data) {
            msg="Could not load image file: "+mapImg;
            textEdit -> append (msg);
        }
        else{
            imshow("Module map", smap );
            int Nrow=smap.rows;
            int Ncol=smap.cols;
            resizeWindow("Module map",Ncol,Nrow);
        }
        QString mapSlopeImg=pathfile+"module_mapSlope.jpeg";
        smapSlope=imread(mapSlopeImg.toStdString(),IMREAD_UNCHANGED);
        if(!smapSlope.data) {
            msg="Could not load image file: "+mapSlopeImg;
            textEdit -> append (msg);
        }
        else{
            imshow("Module Slope-Deviation map", smapSlope );
            int Nrow=smap.rows;
            int Ncol=smap.cols;
            resizeWindow("Module Slope-Deviation map",Ncol,Nrow);
        }
    }
}


//save the measurement file
void VISfieldR::saveMeasFile(){
    QString comment,strg;
    Qt::CheckState state;
    int value;
    double dvalue,Pjp[3];
    fileInfo = QFileDialog::getSaveFileName(this, tr("Choose the File to save"),
                                            pathfile,
                                            tr("infoMeasurement (*.txt)"));
    QFile file(fileInfo);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)){
        msg="Can not open "+fileInfo;
        textEdit -> append (msg);
        return;
    }
    cout<<"saving infoMeasurement.txt"<<"\n";
    QTextStream stream ( &file );
    strg=plant_code -> text();
    stream<<"infoPlant"<<"\t"<<strg<<"\n";
    strg=module_code -> text();
    stream<<"infoModul"<<"\t"<<strg<<"\n";
    strg=FileInfo -> text();
    stream<<"panelType"<<"\t"<<strg<<"\n";
    dvalue=QDSBlengthUP -> value();
    stream<<"moduleLengthUP(mm)"<<"\t"<<dvalue<<"\n";
    dvalue=QDSBlengthDW -> value();
    stream<<"moduleLengthDW(mm)"<<"\t"<<dvalue<<"\n";
    dvalue=dSB_Zc -> value();
    stream<<"distance(mm)"<<"\t"<<dvalue<<"\n";
    stream<<"baseName"<<"\t"<<bNf<<"\n";
    stream<<"imgType"<<"\t"<<eImg<<"\n";
    value=sB_gLup->value();
    stream<<"gLup"<<"\t"<<value<<"\n";
    value=sB_gLdown->value();
    stream<<"gLdw"<<"\t"<<value<<"\n";
    stream<<"frameMIN"<<"\t"<<Nmin<<"\n";
    stream<<"frameMAX"<<"\t"<<Nmax<<"\n";
    dvalue=dSB_mt -> value();
    stream<<"m(deg/sec)"<<"\t"<<dvalue<<"\n";
    dvalue=dSB_q -> value();
    stream<<"q(deg)"<<"\t"<<dvalue<<"\n";
    dvalue=dSB_vjp ->value();
    stream<<"vertex2jp"<<"\t"<<dvalue<<"\n";
    Pjp[0]=doubleSpinBox_xR->value();
    Pjp[1]=doubleSpinBox_yR->value();
    Pjp[2]=doubleSpinBox_zR->value();
    stream<<"JPright"<<"\t"<<Pjp[0]<<"\t"<<Pjp[1]<<"\t"<<Pjp[2]<<"\n";
    Pjp[0]=doubleSpinBox_xL->value();
    Pjp[1]=doubleSpinBox_yL->value();
    Pjp[2]=doubleSpinBox_zL->value();
    stream<<"JPleft"<<"\t"<<Pjp[0]<<"\t"<<Pjp[1]<<"\t"<<Pjp[2]<<"\n";
    if(iFull0Half1==0)
        stream<<"scan"<<"\t"<<"full"<<"\n";
    else
        stream<<"scan"<<"\t"<<"half"<<"\n";
    state = checkBox_fullLength -> checkState();
    if( state == Qt::Checked )
        stream<<"fullLength"<<"\t"<<"yes"<<"\n";
    else
        stream<<"fullLength"<<"\t"<<"no"<<"\n";
    state = checkBox_useOldPoints-> checkState();
    if( state == Qt::Checked )
        stream<<"useOldPoints"<<"\t"<<"yes"<<"\n";
    else
        stream<<"useOldPoints"<<"\t"<<"no"<<"\n";
    state = checkBox_x0here-> checkState();
    if( state == Qt::Checked ){
        stream<<"x0here"<<"\t"<<"yes"<<"\n";
        if(!imgNameVertex.isEmpty())
            stream<<"NimgRef"<<"\t"<<imgNameVertex<<"\n";
    }else{
        stream<<"x0here"<<"\t"<<"no"<<"\n";
        if(NimgRef>0)
            stream<<"NimgRef"<<"\t"<<NimgRef<<"\n";
    }
    if(Np==8 || Np==6){
        for(int i=0; i<Np; i++){
            comment="px["+QString::number(i)+"]";
            stream << comment<<"\t"<< px[i][0]<<"\t"<<px[i][1]<<"\n";
        }
        stream<<"chi2"<<"\t"<<chi2<<"\n";
    }
    file.close();
}


void VISfieldR::updateN(){
    Nmin=sB_Nmin -> value();
    sB_Nimg -> setMinimum(Nmin);
    Nmax=sB_Nmax -> value();
    sB_Nimg -> setMaximum(Nmax);
}


void VISfieldR::viewNmin(){
    sB_Nimg -> setValue(Nmin);
}


void VISfieldR::viewNmax(){
    sB_Nimg -> setValue(Nmax);
}


void VISfieldR::viewImgRef(){
    NimgRef=sB_NimgRef->value();
    sB_Nimg -> setValue(NimgRef);
}



void VISfieldR::viewFrameX(){
    if(ioBrek==1)
        return;
    double DthetaAbs,thetaAbs;
    Nimg=sB_Nimg -> value();
    printf("************* viewFrameX() Nimg= %d\n",Nimg);
    QString imgName;
    imgName=bNf+QString::number(Nimg).rightJustified(4, '0')+eImg;
    imgName=pathfile+imgName;
    lineEdit_img ->setText(imgName.section('/', -4));
    int iMeth=comboBox_method->currentIndex();
    if(iMeth==0){
        doubleSpinBox_offset->setEnabled(true);
        double stepTheta=dSB_thetaStep->value();
        double offset=doubleSpinBox_offset->value();
        DthetaAbs=(NimgRef-Nimg)*stepTheta+offset;
        thetaAbs=90.+DthetaAbs;
    }
    else{
        doubleSpinBox_offset->setEnabled(false);
        valT=JPGtime(imgName);
        double mt=  dSB_mt -> value();
        double q=  dSB_q -> value();
        //printf("JPGtime: mt=%f q=%f\n",mt,q);
        thetaAbs=mt*valT+q;
        if(thetaAbs>0.)
            DthetaAbs=thetaAbs-90.;
        else
            DthetaAbs=thetaAbs+90.;
    }
    lineEdit_thetaAbs -> setText(QString::number(thetaAbs));
    printf(" DthetaAbs(deg)=%f\n",DthetaAbs);
    lineEdit_theta -> setText(QString::number(DthetaAbs));
    thetaImg=DthetaAbs;
    thetaImg=thetaImg*pig/180.;
    showframe(imgName);
}


double VISfieldR::JPGtime(QString imgName){
    int valH,valM,valS,valCS;
    double valTime;
    QString command,Vline;
    //read date and time of JPG image and write the values in info.txt
    command="exiftool -s3 -DateTimeOriginal -SubSecTime "+imgName+" > "+pathbase+"info.txt";
    system(command.toStdString().c_str());
    fflush(stdout);
    waitKey(10);
    QFile file(pathbase+"info.txt");
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        cout << "ERROR opening: "<<(pathbase+"info.txt").toStdString() <<"\n";
    QTextStream in(&file);
    Vline = in.readLine();
    valH=(Vline.midRef(11,2)).toInt();
    valM=(Vline.midRef(14,2)).toInt();
    valS=(Vline.midRef(17,2)).toInt();
    Vline = in.readLine();
    valCS=Vline.toInt();
    valTime=double(valH*3600)+double(valM*60)+double(valS)+double(valCS)/100.;
    printf("h=%d m=%d s=%d.%d T(sec)=%f\n",valH,valM,valS,valCS,valTime);
    lineEdit_t1 -> setText(QString::number(valH)+":"+QString::number(valM)+":"+QString::number(double(valS)+double(valCS)/100.,'f',3));
    file.close();
    return(valTime);
}


void VISfieldR::showframe(QString spath){
    QString filetmp;
    printf("->showFrame() with ...\n");
    cout << "\timg= "<<spath.toStdString()<<"\n";
    cout << "\tthetaImg= "<<thetaImg<<"\n";
    cout << "\tNp= "<<Np<<"\n";
    cout << "\tNimg= "<<Nimg<<"\n";
    double v2jp=dSB_vjp ->value();
    double sundiv=DSBsundiv -> value();
    sundiv=sundiv/1000.0;
    double sundecli=DSBsun_decli -> value();//deg
    radius=QDSBraggio -> value();
    double sundecliMAX=acos(tan(sundiv)*sqrt((focal/radius)*(focal/radius)-1))*180.0/pig;
    double tg2p=sundiv/cos(sundecli/180.0*pig);
    tg2p=tg2p*tg2p;
    if(sundecli>sundecliMAX){
        msg="because of the large value of Sun declination\nnot all solar radiation can be captured!";
        QMessageBox msgBox;
        msgBox.setText(msg);
        msgBox.exec();
    }
    imgDisplayed=imread(spath.toStdString(),IMREAD_UNCHANGED);
    if(!imgDisplayed.data) {
        msg="Could not load image file: "+spath;
        textEdit -> append (msg);
        return;
    }
    //undistort(imgDisplayed,img,cameraMatrix,distCoeffs);
    img=imgDisplayed.clone();
    cvtColor(img, img8b, COLOR_RGB2GRAY);
    filetmp=pathbase+"img.JPG";
    imwrite(filetmp.toStdString().c_str(),img8b);
    imgDisplayed=img.clone();
    Width    = img.cols;
    Height   = img.rows;
    lineEdit_width  ->setText(QString::number(Width));
    lineEdit_height ->setText(QString::number(Height));
    sB_j->setMaximum(Width);
    sB_i->setMaximum(Height);
    if (img.depth() == CV_8U){
        sB_gLup -> setMaximum(255);
        sB_gLdown -> setMaximum(255);
        lineEdit_depth -> setText("8 bit");
    }else if (img.depth() == CV_16U){
        sB_gLup -> setMaximum(65535);
        sB_gLdown -> setMaximum(65535);
        lineEdit_depth -> setText("16 bit");
    }
    lineEdit_channels -> setText(QString::number(img.channels()));
    //check gray threshold
    Qt::CheckState state1;
    state1 = checkBoxGray -> checkState();
    if(iFull0Half1==1)
        soglia=sB_gLdown -> value();
    else{
        if(Nimg<(Nmin+Nmax)/2)
            soglia=sB_gLup -> value();
        else
            soglia=sB_gLdown -> value();
    }
    if( state1 == Qt::Checked ) {
        for(int k=0;k<Height;k++){
            for(int l=0;l<Width;l++){
                if(img8b.at<uchar>(k,l)<=soglia) {
                    imgDisplayed.at<Vec3b>(k,l)[0]=255;
                    imgDisplayed.at<Vec3b>(k,l)[1]=0;
                    imgDisplayed.at<Vec3b>(k,l)[2]=0;
                }
            }
        }
    }
    if(Np>=2)//traccia la linea tra i due target joint-pivot
        line(imgDisplayed,Point(int(px[0][0]+cx),int(px[0][1]+cy)),Point(int(px[1][0]+cx),int(px[1][1]+cy)),Scalar(255,0,0), 3);
    if(iSetRef==1){
        //setGap();
        liner(1,thetaImg,0.,-length/2.,0.,0.,length/2.,0.);//asse rotazione
        liner(2,thetaImg,0.,-length/2.,focal-v2jp,0.,length/2.,focal-v2jp);//asse focale
        //camera position in the turned parabola reference frame
        Mod2Dro(0.,-thetaImg,0.,Pc[0],Pc[1],Pc[2]);//serve il segno -
        xCa=Vservice[0];
        yCa=Vservice[1];
        zCa=Vservice[2]+v2jp;
        xP=Vservice[0];
        lineEdit_xCam->setText(QString::number(xP,'f', 1 ));
        zP=0.25/focal*xP*xP-v2jp;
        dview=sqrt((xP-Pc[0])*(xP-Pc[0])+(zP-Pc[2])*(zP-Pc[2]));
        lineEdit_dview->setText(QString::number(int(dview)));
        imgDimension(tg2p);//limits of HCE and solar spot -> Magnification
        lineEdit_solSpot->setText(QString::number(int(fabs(xsMax-xsMin))));
        if(sundecli>sundecliMAX){
            xrMin=xsMin;
            xrMax=xsMax;
        }
        //traccia linee immagine e spot solare su immagine aerea dato x
        if(xrMin<-xmax) xrMin=-xmax;
        if(xrMin>xmax) xrMin=xmax;
        if(xrMax<-xmax) xrMax=-xmax;
        if(xrMax>xmax) xrMax=xmax;
        if(xsMin<-xmax) xsMin=-xmax;
        if(xsMin>xmax) xsMin=xmax;
        if(xsMax<-xmax) xsMax=-xmax;
        if(xsMax>xmax) xsMax=xmax;
        //traccia bordi modulo
        zP=0.25/focal*xmax*xmax-v2jp;
        liner(1,thetaImg,xmax,-lengthDW/2.,zP,xmax,lengthDW/2.,zP);//+xmax
        liner(1,thetaImg,-xmax,-lengthUP/2.,zP,-xmax,lengthUP/2.,zP);//-xmax
        //traccia separazione pannelli
        if(iPanel==1){
            double xmedio=QDSBinnerXmax->value();
            xmedio=xmedio+QDSBouterXmin->value();
            xmedio=xmedio/2.;
            zP=0.25/focal*xmedio*xmedio-v2jp;
            liner(1,thetaImg,xmedio,-lengthDW/2.,zP,xmedio,lengthDW/2.,zP);//+xmedio
            liner(1,thetaImg,-xmedio,-lengthUP/2.,zP,-xmedio,lengthUP/2.,zP);//-xmedio
        }
        //traccia linea ricevitore xrMIN
        zP=0.25/focal*xrMin*xrMin-v2jp;
        liner(3,thetaImg,xrMin,-length/2.,zP,xrMin,length/2.,zP);//xrMin
        //traccia linea ricevitore xrMAX
        zP=0.25/focal*xrMax*xrMax-v2jp;
        liner(3,thetaImg,xrMax,-length/2.,zP,xrMax,length/2.,zP);//xrMax
        // traccia linea spot solare xsMIN
        zP=0.25/focal*xsMin*xsMin-v2jp;
        liner(3,thetaImg,xsMin,-length/2.,zP,xsMin,length/2.,zP);//xsMin
        // traccia linea spot solare xsMAX
        zP=0.25/focal*xsMax*xsMax-v2jp;
        liner(3,thetaImg,xsMax,-length/2.,zP,xsMax,length/2.,zP);//xsMax
        //plot reference points
        for(int i=0;i<Np;i++){
            pointer(1,thetaImg,P[i][0],P[i][1],P[i][2]);
        }
        //analisi nella fascia di ricerca
        int intensity;
        int irs=Nimg-Nmin;
//        int irs=int((xP+xmax)/xStep);
//        if(iFull0Half1==1)
//            irs=int((xP)/xStep);
//        printf("xP=%f xStep=%f irs=%d\n",xP,xStep,irs);
        Nirs=max(Nirs,irs);
        int jrs;
        double icont,isum,Dpx,ncont,nifOK,s1,s2;
        double xSearchMin=max(xrMin-(xrMax-xrMin)/2.,-xmax);
        double xSearchMax=min(xrMax+(xrMax-xrMin)/2.,xmax);
        printf("xrMin=%f xrMax=%f\nxSearchMin=%f xSearchMax=%f\n",xrMin,xrMax,xSearchMin,xSearchMax);
        ROI(thetaImg,xSearchMin,xSearchMax);//fascia di ricerca per valutazione baricentro
        double mjte=mjt;
        double cjte=cjt;
        double mjbe=mjb;
        double cjbe=cjb;
        double mile=mil;
        double cile=cil;
        double mire=mir;
        double cire=cir;
        double intFac=0.;
        double pxroi,delta;
        double fatCorr;
        s1=0.;
        s2=0.;
        ROI(thetaImg,xsMin,xsMax);//fascia spot solare
        for(int j=0;j<NpxX;j++){
            icont=0.;
            isum=0.;
            ncont=0.;
            nifOK=0.;
            pxroi=(j*mjt+cjt+j*mjb+cjb)/2.;
            fatCorr=(jRight*mjb+cjb-jRight*mjt-cjt)/(j*mjb+cjb-j*mjt-cjt);//correction for oblique view
            Dpx=sqrt((j-jLeft)*(j-jLeft)+(pxroi-iLeft)*(pxroi-iLeft));
            //printf("fatCorr=%f jrs=%d jrs_new=%d\n",fatCorr,int(Dpx/Lpx*double(Njrs)),int(Dpx/Lpx*double(Njrs)*fatCorr));
            jrs=int(Dpx/Lpx*double(Njrs)*fatCorr);
//            if(j==0 || j==NpxX-1)
//                printf("jrs@j=%d = %d\n",j,jrs);
            for(int i=0;i<NpxY;i++){
                if(i > int(j*mjte+cjte+0.5) &&
                        i < int(j*mjbe+cjbe+0.5) &&
                        j > int(i*mile+cile+0.5) &&
                        j < int(i*mire+cire+0.5)   ){
                    intensity= img8b.at<uchar>(i, j);
//                    if(j==3700)
//                        printf("i=%d intensity=%d\n",i,intensity);
                    if(intensity <= soglia){
                        icont++;
                        isum=isum+double(i);
                    }
                }
                if(i > int(j*mjt+cjt+0.5) &&
                        i < int(j*mjb+cjb+0.5) &&
                        j > int(i*mil+cil+0.5) &&
                        j < int(i*mir+cir+0.5)   ){
                    ncont++;
                    s1++;
                    intensity= img8b.at<uchar>(i, j);
                    if(intensity <= soglia){
                        nifOK++;
                        s2++;
                    }
                }
            }
            //printf("icont=%f ncont=%f\n",icont,ncont);
            //waitKey(0);
            if(icont>1. && ncont>1.){
                delta=dpximg*(isum/icont-pxroi)-((xrMax+xrMin)/2.-xP);
                mappa.at<double>(irs,jrs,0)++;
                mappa.at<double>(irs,jrs,2)=mappa.at<double>(irs,jrs,2)+delta/2.0/Magne/focal;//arctangent deviation
                mappa.at<double>(irs,jrs,3)=mappa.at<double>(irs,jrs,3)+nifOK/ncont;
                mappa.at<double>(irs,jrs,4)=cos(thetaImg);//peso della riga i-esima di mappa
                mappa.at<double>(irs,jrs,5)=dpximg/2.0/Magne/focal;//sensitivity
            }
        }
        //printf("s1=%f s2=%f\n",s1,s2);
        intFac=s2/s1;
        lineEdit_intFat->setText(QString::number(100.*intFac,'f', 1 ));
        if(xP>=-xmax && xP<=xmax){
            IntFact[Nimg-Nmin][0]=1.;
            IntFact[Nimg-Nmin][1]=intFac;
        }
        else
            IntFact[Nimg-Nmin][0]=0.;
        int imgProcessed=0;
        double sum=0.;
        for(int i=0;i<=Nmax-Nmin;i++){
            if(IntFact[i][0]>0.){
                sum=sum+IntFact[i][1];
                imgProcessed++;
            }
        }
        if(imgProcessed>0)
            lineEdit_IFmean->setText(QString::number(100.*sum/double(imgProcessed),'f', 1 ));
        else
            lineEdit_IFmean->setText(QString::number(0.,'f', 1 ));
        QFile ftP(fileResult);//per salvataggio risultati
        int iOKopen=1;
        if(!ftP.open(QIODevice::ReadWrite | QIODevice::Text))
            iOKopen=0;
        QTextStream stream ( &ftP );
        if(iOKopen==1){
            stream<<Nimg<<"\t"<<xP<<"\t"<<thetaImg<<"\t"<<intFac<<"\n";
            ftP.close();
        }
    }
    imshow("Image", imgDisplayed );
    filetmp=pathbase+"imgDisplayed.JPG";
    imwrite(filetmp.toStdString().c_str(),imgDisplayed);
    if(cropInfo[0]==1){
        Rect myROI(cropInfo[1],cropInfo[2],cropInfo[3],cropInfo[4]);//ROI
        Cimg = imgDisplayed(myROI);// cropped image
        imshow( "ROI",Cimg);
    }
}


void VISfieldR::process(){
    //inzialization
    Nirs=0;
    double panel_x[2*Nraw],PerSurf[7][20][10],weightRow[Nraw];
    if(iFull0Half1==0){
        int iC_Phang=Nraw;
        panel_x[iC_Phang]= xmin;
        panel_x[iC_Phang-1]= -xmin;
        panel_x[iC_Phang+1]= QDSBinnerXmax -> value();
        panel_x[iC_Phang-2]=-panel_x[iC_Phang+1];
        if(iPanel == 1){
            panel_x[iC_Phang+2]= QDSBouterXmin -> value();
            panel_x[iC_Phang+3]= QDSBouterXmax -> value();
            panel_x[iC_Phang-3]=-panel_x[iC_Phang+2];
            panel_x[iC_Phang-4]=-panel_x[iC_Phang+3];
        }
        else if(iPanel == 2){
            panel_x[iC_Phang+2]= QDSBmiddleXmin -> value();
            panel_x[iC_Phang+3]= QDSBmiddleXmax -> value();
            panel_x[iC_Phang-3]=-panel_x[iC_Phang+2];
            panel_x[iC_Phang-4]=-panel_x[iC_Phang+3];
            panel_x[iC_Phang+4]= QDSBouterXmin -> value();
            panel_x[iC_Phang+5]= QDSBouterXmax -> value();
            panel_x[iC_Phang-5]=-panel_x[iC_Phang+4];
            panel_x[iC_Phang-6]=-panel_x[iC_Phang+5];
        }
        for(int i=0;i<2*Nraw;i++){
           printf("panel_x[%d] = %f\n",i,panel_x[i]);
        }
        //Row weights are given by the relative cross-section
        if(iPanel == 0){
            weightRow[0]=0.5;
            weightRow[1]=0.5;
        }
        else if(iPanel ==1){
            weightRow[0]=(panel_x[7]-panel_x[6])/(panel_x[5]-panel_x[4]+panel_x[7]-panel_x[6])/2.;//outer
            weightRow[1]=(panel_x[5]-panel_x[4])/(panel_x[5]-panel_x[4]+panel_x[7]-panel_x[6])/2.;//inner (of the semiparabola)
            weightRow[2]=weightRow[1];//inner
            weightRow[3]=weightRow[0];//outer
        }
        else if(iPanel ==2){
            weightRow[0]=(panel_x[11]-panel_x[10])/(panel_x[11]-panel_x[10]+panel_x[9]-panel_x[8]+panel_x[7]-panel_x[6])/2.;//outer
            weightRow[1]=(panel_x[9]-panel_x[8])  /(panel_x[11]-panel_x[10]+panel_x[9]-panel_x[8]+panel_x[7]-panel_x[6])/2.;//middle
            weightRow[2]=(panel_x[7]-panel_x[6])  /(panel_x[11]-panel_x[10]+panel_x[9]-panel_x[8]+panel_x[7]-panel_x[6])/2.;//inner
            weightRow[3]=weightRow[2];//inner
            weightRow[4]=weightRow[1];//middle
            weightRow[5]=weightRow[0];//outer
        }
    }
    else{
        panel_x[0]= xmin;
        panel_x[1]= QDSBinnerXmax -> value();
        if(iPanel == 1){
            panel_x[2]= QDSBouterXmin -> value();
            panel_x[3]= QDSBouterXmax -> value();
        }
        else if(iPanel == 2){
            panel_x[2]= QDSBmiddleXmin -> value();
            panel_x[3]= QDSBmiddleXmax -> value();
            panel_x[4]= QDSBouterXmin -> value();
            panel_x[5]= QDSBouterXmax -> value();
        }
        for(int i=0;i<6;i++)
            panel_x[i]=panel_x[i];
        if(iPanel == 0){
            weightRow[0]=1.0;
            weightRow[1]=0.0;
        }
        else if(iPanel ==1){
            weightRow[0]=(panel_x[1]-panel_x[0])/(panel_x[1]-panel_x[0]+panel_x[3]-panel_x[2]);
            weightRow[1]=(panel_x[3]-panel_x[2])/(panel_x[1]-panel_x[0]+panel_x[3]-panel_x[2]);
        }
        else if(iPanel ==2){
            weightRow[0]=(panel_x[1]-panel_x[0])/(panel_x[1]-panel_x[0]+panel_x[3]-panel_x[2]+panel_x[5]-panel_x[4]);
            weightRow[1]=(panel_x[3]-panel_x[2])/(panel_x[1]-panel_x[0]+panel_x[3]-panel_x[2]+panel_x[5]-panel_x[4]);
            weightRow[2]=(panel_x[5]-panel_x[4])/(panel_x[1]-panel_x[0]+panel_x[3]-panel_x[2]+panel_x[5]-panel_x[4]);
        }
    }
    for(int i=0;i<Nraw;i++){
       printf("weightRow[%d] = %f\n",i,weightRow[i]);
    }
    //reset mappa
    for(int jrs=0;jrs<NjMappa;jrs++){
        for(int irs=0;irs<NiMappa;irs++){
            for(int k=0;k<4;k++)
                mappa.at<double>(irs,jrs,k)=0.;
        }
    }
    //open file where write xP
    QFile fileXP(pathfile+"xP.txt");
    if (!fileXP.open(QIODevice::WriteOnly | QIODevice::Text)){
        msg="Can not open "+pathfile+"xP.txt";
        textEdit -> append (msg);
        return;
    }
    QTextStream streamXP(&fileXP);
    irsStep=1;
    double xPold;
    //process loop
    for(int nf=Nmin;nf<=Nmax;nf++){
        //printf("Nframe= %d\t",nf);
        sB_Nimg -> setValue(nf);//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< view and analyze frames
        streamXP << nf <<"\t" << thetaImg <<"\t" << xCa <<"\t" << yCa<<"\t" << zCa << "\n";
        if(nf==Nmin)
            xPold=xP;
        if(nf==Nmin+1){
            if(xP<xPold)
                irsStep=-1;
        }
        waitKey(10);
    }
    fileXP.close();
    //end processing
    printf("********** END PROCESSING *********\n");

    //normalization
    for(int irs=0;irs<=Nirs;irs++){
        for(int jrs=0;jrs<=Njrs;jrs++){
            if(mappa.at<double>(irs,jrs,0)>1.){
                mappa.at<double>(irs,jrs,1)=55.;
                mappa.at<double>(irs,jrs,2)=mappa.at<double>(irs,jrs,2)/mappa.at<double>(irs,jrs,0);
                mappa.at<double>(irs,jrs,3)=mappa.at<double>(irs,jrs,3)/mappa.at<double>(irs,jrs,0);
                if(mappa.at<double>(irs,jrs,3) > 0.999)
                    mappa.at<double>(irs,jrs,1)=255.;
            }
        }
    }
    iProcess=1;

    //plot map
    plotMap();

    //******** data analysis for each panel composing the collector
    char command[200];
    int Npanel=int(length/(panel_width)+0.5);//N panel along one raw
    int width_true=Njrs;
    int nheight=Nirs;
    int ipanel,iraw,iHangDW=0,iHangUP=0,iHang0=0;
    double sumx,sumxx,sumxy,sumy,sumyy,sumErr,sumIF,N,a,b,sampled,aligned,mu,rms,
            y_left,y_right,ErrAverage,IntFact,IntFactTot,suax,suaxx,suay,suaxy,Nalpha;
    double unbelt=50.;//mm
    double Eradius=62.5;//mm
    double sundiv=0.00473;//rad
    double freccia=0.;//mm
    double sundecli=0.;//deg
    double gap=gapDW;
    if(gap < 0.0){
        msg="Error!\ngap-value has to be POSITIVE";
        QMessageBox msgBox;
        msgBox.setText(msg);
        msgBox.exec();
    }
    dpximg=lengthDW/double(Njrs);
    int Npoint,NpOK,j1,j2,jc,i1,i2,ierr;
    int j_step=int((panel_width+gap)/dpximg+0.5);//step in pixel between one panel and the following
    int j_stepD=int((panel_width/2.0-unbelt)/dpximg+0.5);//panel halfwidth to analyse
    int j_offeset=0;
    int ipc_lim=int(double(Npanel)/2.0);
    if((Npanel-2*int(double(Npanel)/2.0+0.5)) == 0){
        j_offeset=int((panel_width+gap)/2.0/dpximg+0.5);
        ipc_lim=int(double(Npanel)/2.0)-1;
    }
    if(iFull0Half1==1)
        msg="      the half-module is composed by "+ QString::number(Npanel*Nraw) +
            " panel disposed on " + QString::number(Nraw) +" raw";
    else
        msg="      the module is composed by "+ QString::number(Npanel*Nraw) +
                " panel disposed on " + QString::number(Nraw) +" raw";
    textEdit -> append (msg);
    msg="        each panel is hang on " + QString::number(N_Phang) + " points";
    textEdit -> append (msg);
    cout<<msg.toStdString()<<"\n";

    IntFactTot=0.0;
    for(int iz=0;iz<20;iz++){
        for(int ix=0;ix<10;ix++)
            corrPhold[iz][ix]=0.0;
    }
    suax=0.0;
    suaxx=0.0;
    suaxy=0.0;
    suay=0.0;
    Nalpha=0.0;
    // if(Nraw>=2){
    //     int im=int((panel_x[1]-panel_x[0])/dpximg+0.5);
    //     line(smap, Point(0,im), Point(width_true-1,im), cvScalar(0,255,0),1);
    //     im=int((panel_x[2]-panel_x[0])/dpximg+0.5);
    //     line(smap, Point(0,im), Point(width_true-1,im), cvScalar(0,255,0),1);
    // }
    // if(Nraw==3){
    //     int im=int((panel_x[3]-panel_x[0])/dpximg+0.5);
    //     line(smap, Point(0,im), Point(width_true-1,im), cvScalar(0,255,0),1);
    //     im=int((panel_x[4]-panel_x[0])/dpximg+0.5);
    //     line(smap, Point(0,im), Point(width_true-1,im), cvScalar(0,255,0),1);
    // }
    // for(int iSide=0;iSide<2*Nraw;iSide++){// draw the panel edges on smap
    //     int Di=int((panel_x[iSide]-panel_x[0])/dpximg+0.5);
    //     Di=nheight/2+segno(double(iSide))*Di;
    //     line(smap, Point(0,Di), Point(width_true-1,Di), Scalar(0,255,0),1);
    // }
     printf("Nraw=%d N_Phang=%d\n",Nraw,N_Phang);
    for(iraw=0;iraw<Nraw;iraw++){
        i1=int((panel_x[0+2*iraw]-panel_x[0])/dpximg+0.5);
        i2=int((panel_x[1+2*iraw]-panel_x[0])/dpximg+0.5);
        if(i2>=nheight) i2=nheight-1;
        iHangUP=iraw*N_Phang;//index of upper hanging points
        if(N_Phang==4){//2
            iHangDW=iHangUP+1;
            iHang0=iHangDW;
        }else if(N_Phang==6){//3
            iHang0=iHangUP+1;
            iHangDW=iHangUP+2;
        }else{
            msg="The current Number of hanging points per panel is not supported!\n"
                "The software has to be implemented for the new number!";
            QMessageBox msgBox;
            msgBox.setText(msg);
            msgBox.exec();
        }
        //printf("iraw= %d\ti1= %d\ti2= %d\tiHangUP= %d\tiHang0= %d\tiHangDW= %d\n",iraw,i1,i2,iHangUP,iHang0,iHangDW);
        for(int ipc=0;ipc<Npanel;ipc++){
            sumx=0.0;
            sumxx=0.0;
            sumxy=0.0;
            sumy=0.0;
            sumyy=0.0;
            sumErr=0.0;
            sumIF=0.0;
            Npoint=0;
            NpOK=0;
            if(i2>=nheight) i2=nheight-1;
            if(ipc <= ipc_lim){
                ipanel=int(double(Npanel)/2.0)+ipc;
                jc=width_true/2+ipc*j_step+j_offeset;
            }
            else{
                ipanel=Npanel-1 - ipc;
                jc=width_true/2-(ipc-int(double(Npanel)/2.0))*j_step-j_offeset;
            }
            j1=jc-j_stepD;
            j2=jc+j_stepD;
            if(j1 < 0) j1=0;
            if(j2 > width_true-1)
                j2=width_true-1;
            printf("ipc=%d ipanel=%d jc=%d j1=%d j2=%d \n",ipc,ipanel,jc,j1,j2);
            line(smap, Point(j1,0), Point(j1,nheight-1), Scalar(0,255,0),1);
            line(smap, Point(j2,0), Point(j2,nheight-1), Scalar(0,255,0),1);
            int jmedio=int(double(j2+j1)/2.0+0.5);
            for(int i=i1;i<i2;i++) {
                for(int j=j1;j<j2;j++) {
                    if(mappa.at<double>(i,j,1)>50.0){
                        if(mappa.at<double>(i,j,1)>250.0) NpOK++;
                        sumx=sumx+(j-jmedio)*dpximg*mappa.at<double>(i,j,4);
                        sumxx=sumxx+(j-jmedio)*(j-jmedio)*dpximg*dpximg*mappa.at<double>(i,j,4);
                        sumxy=sumxy+mappa.at<double>(i,j,2)*(j-jmedio)*dpximg*mappa.at<double>(i,j,4);
                        sumy=sumy+mappa.at<double>(i,j,2)*mappa.at<double>(i,j,4);
                        sumyy=sumyy+mappa.at<double>(i,j,2)*mappa.at<double>(i,j,2)*mappa.at<double>(i,j,4);
                        sumErr=sumErr+mappa.at<double>(i,j,5)*mappa.at<double>(i,j,4);
                        sumIF=sumIF+mappa.at<double>(i,j,3)*mappa.at<double>(i,j,4);
                        Npoint++;
                    }
                }
            }
            N=double(Npoint);
            mu=sumy/N;
            rms=sqrt(sumyy/N);
            ErrAverage=sumErr/N;
            a = (N*sumxy-sumx*sumy)/(N*sumxx-sumx*sumx);
            b = sumy/N - a*sumx/N;
            //      SEa = ErrAverage/sqrt(sumxx-sumx*sumx/N);
            //      SEb = ErrAverage*sqrt(1.0/N+sumx/N*sumx/N/(sumxx-sumx*sumx/N));
            sampled=N/double((i2-i1)*(j2-j1));
            aligned=double(NpOK)/N;
            ErrAverage=sumErr/N;
            IntFact=sumIF/N;
            IntFactTot=IntFactTot+IntFact*weightRow[iraw];
                 printf("***** Raw %d Panel N. %d average=%e rms=%e <Err>=%e\n",iraw+1,ipanel+1,mu,rms,ErrAverage);
                 printf("j1=%d 12=%d\n",j1,j2);
                 printf("a = %e   b = %e \n",a,b);
                 printf(" Sampled-surface fraction: %4.3f  ... of which rigthly aligned: %4.3f\n",sampled,aligned);
            PerSurf[iraw][ipanel][0]=sampled;
            PerSurf[iraw][ipanel][1]=aligned;
            PerSurf[iraw][ipanel][2]=mu;
            PerSurf[iraw][ipanel][3]=rms;
            PerSurf[iraw][ipanel][4]=ErrAverage;
            PerSurf[iraw][ipanel][5]=IntFact;
            y_left =-(Phold[iHangDW][0]-Phold[iHang0][0])*(b-a*Phold[iHangDW][1]);//suggested corrections
            y_right=-(Phold[iHangDW][0]-Phold[iHang0][0])*(b+a*Phold[iHangDW][1]);
            corrPhold[ipanel*2  ][iHangDW]=y_left;
            corrPhold[ipanel*2+1][iHangDW]=y_right;
            if(N_Phang==6){//3
                y_left =-(Phold[iHangUP][0]-Phold[iHang0][0])*(b-a*Phold[iHangUP][1]);//suggested corrections
                y_right=-(Phold[iHangUP][0]-Phold[iHang0][0])*(b+a*Phold[iHangUP][1]);
                corrPhold[ipanel*2  ][iHangUP]=y_left;
                corrPhold[ipanel*2+1][iHangUP]=y_right;
            }
            double x1=Phold[iHang0][1]+(ipanel+0.5)*(panel_width+gap);
            double y1=-(b-a*Phold[iHang0][1]);
            double x2=Phold[iHang0][1]+(ipanel+0.5)*(panel_width+gap);
            double y2=-(b+a*Phold[iHang0][1]);
            matAlpha[ipanel*2  ][iraw]=y1;
            matAlpha[ipanel*2+1][iraw]=y2;
            //if(ipanel>0 && ipanel<Npanel-1){
                Nalpha=Nalpha+2.0;
                suax=suax+x1+x2;
                suaxx=suaxx+x1*x1+x2*x2;
                suaxy=suaxy+x1*y1+x2*y2;
                suay=suay+y1+y2;
            //}
        }
    }
    printf("suax=%e suaxx=%e suaxy=%e suay=%e\n",suax,suaxx,suaxy,suay);

    if(N_Phang==6){//3 adjust corrections to be always >=0
        double Delta;
        for(int iraw=0;iraw<Nraw;iraw++){
            iHangUP=iraw*N_Phang;//index of upper hanging points
            iHang0=iHangUP+1;
            iHangDW=iHangUP+2;
            for(int iz=0;iz<Npanel*2;iz++){
                if(corrPhold[iz][iHangUP]<0)
                    Delta=fabs(corrPhold[iz][iHangUP]);
                else
                    Delta=fabs(corrPhold[iz][iHangDW]);
                corrPhold[iz][iHangUP]=corrPhold[iz][iHangUP]+Delta;
                corrPhold[iz][iHang0 ]=corrPhold[iz][iHang0 ]+Delta;
                corrPhold[iz][iHangDW]=corrPhold[iz][iHangDW]+Delta;
            }
        }
    }

    double alpha_medio=suay/Nalpha;
    a = (Nalpha*suaxy-suax*suay)/(Nalpha*suaxx-suax*suax);
    b = suay/Nalpha - a*suax/Nalpha;
    double torsion=a*length;
    IntFactTot=IntFactTot/double(Npanel);

    msg="         px dimension of the module_map = " + QString::number(dpximg,'f',1) +" mm";
    textEdit -> append(msg);
    msg="          module mean-intercept_factor = " + QString::number(IntFactTot,'f',3);
    textEdit -> append(msg);
    msg="           alpha_mean (mrad) = " + QString::number(alpha_medio*1000.0,'f',3) +
            "\n           module torsion (mrad) = "+ QString::number(torsion*1000.0,'f',3);
    textEdit -> append(msg);

    //correction dat file<<<<<<<<<<<<<<<<<<<
     QString filesave=pathbase+"corrections.dat";
     ofstream myfile;
    myfile.open (filesave.toStdString().c_str());
    if(iFull0Half1==1)
        myfile<<"#z(mm) outR1 midR1 innR1 outR2 midR2 innR2"<<"\n";
    else
        myfile<<"#z(mm) outR1 midR1 innR1 outR2 midR2 innR2 innR3 midR3 outR3 innR4 midR4 outR4"<<"\n";
     int mp=-1;
    for(int iz=0;iz<Npanel*2;iz++){
        myfile << panel_width/2.+mp*Phold[iHang0][1]+(int(iz/2))*(panel_width+gap)  << "\t";
        mp=-mp;
        for(int ix=0;ix<N_Phang*Nraw;ix++)
            myfile << corrPhold[iz][ix] << "\t";
        myfile << "\n";
    }
    myfile.close();

    //alpha dat file<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    filesave=pathbase+"alpha.dat";
    ofstream myfilea;
    myfilea.open (filesave.toStdString().c_str());
    myfilea<<"#z(mm)\talphaR1\talphaR2\talphaR3\talphaR4"<<"\n";
    mp=-1;
    for(int iz=0;iz<Npanel*2;iz++){
        myfilea << panel_width/2.+mp*Phold[iHang0][1]+(int(iz/2))*(panel_width+gap)  << "\t";
        mp=-mp;
        for(int ix=0;ix<Nraw;ix++)
            myfilea << matAlpha[iz][ix]*1000. << "\t";
        myfilea << "\n";
    }
    myfilea.close();

    //IF dat file<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    filesave=pathbase;
    filesave+="IntFact.dat";
    myfile.open (filesave.toStdString().c_str());
    for(ipanel=0;ipanel<Npanel;ipanel++){
        sumx=0.;
        myfile << ipanel +1;
        for(iraw=0;iraw<Nraw;iraw++){
            sumx=sumx+PerSurf[iraw][ipanel][5];
            myfile <<"\t"<< PerSurf[iraw][ipanel][5] ;
        }
        myfile << "\t"<<sumx/double(Nraw)<<"\n";
    }
    myfile << "#InterceptFactor_mean " << IntFactTot << "\n";
    myfile.close();

    // copy file.dat in pathbase
    pippo = sprintf(command,"cp %s*.dat %s/ ",pathbase.toStdString().c_str(),pathfile.toStdString().c_str());
    ierr=system(command);

    //***** PLOTS !  with gnuplot **************************
    //correction plot
    if(iPanel==0)
        GnuPlot("plot1a.sh");
    else if(iPanel==1){
        if(iFull0Half1==1)
            GnuPlot("plot1b.sh");
        else
            GnuPlot("plot1bFull.sh");
    }
    else if(iPanel==2)
        GnuPlot("plot1c.sh");

    //alpha plot
    if(iFull0Half1==1)
        GnuPlot("plot2.sh");
    else
        GnuPlot("plot2Full.sh");
    //IF map plot
    if(iPanel==0)
        GnuPlot("plot3a.sh");
    else if(iPanel==1){
        if(iFull0Half1==1)
            GnuPlot("plot3b.sh");
        else
            GnuPlot("plot3bFull.sh");
    }
    else if(iPanel==2)
        GnuPlot("plot3c.sh ");
    // ***** end PLOT *************************************

    // copy file.png in sample folder
    pippo = sprintf(command,"cp %s*.png %s",pathbase.toStdString().c_str(),pathfile.toStdString().c_str());
    ierr=system(command);

//    //convert ps in jpg
//    pippo = sprintf(command,"%s/ps2jpg.sh ",pathbase.toStdString().c_str());
//    ierr=system(command);
//    if(ierr!=0) {
//        msg="not found "+pathbase + "/ps2jpg.sh"+
//                "\n ps2ipg image not executed!!!";
//        QMessageBox::critical(0, "Error", msg);
//    }

    //creating the CSV file for Tellico database
    QString fileInfo = FileInfo -> text();
    QString plant=plant_code ->text();
    QString module=module_code ->text();
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    string v=",";
    QFileInfo info1(pathfile);
    QDateTime qdt = info1.birthTime();
    QString qdtst=qdt.toString( "yyyy-MM-dd") ;
    filesave=pathbase;
    filesave+="data.csv";
    myfile.open (filesave.toStdString().c_str());
    myfile << "ID,Plant,Module,Lmodule,PanelType,Lpanel,gap,Focal,r_rec,r_env,Sdiv,Sdec,Ttr,P_hang,"
           << "image,pxMap,IFmap,v-c,POI_left,POI_right,grayTreshold,UnBelt,IF_map,"
           << "IF_plot,Corrections_plot,AngCorr_plot,Date Created,Date Modified\n";
    myfile<<"0"<<v<<plant.toStdString()<<v<< module.toStdString()<<v<<length<<v<<PanelType.toStdString()<<
            v<<panel_width<<v<<gap<<v<<focal<<v<<radius<<v<<Eradius<<v<<sundiv<<v<<sundecli<<v<<freccia<<
            v<<fileInfo.toStdString()<<v<<Width<<"x"<<Height<<v<<dpximg<<v<<width_true<<"x"<<nheight<<v<<dview<<
            v<<"NA"<<v<<"NA"<<v<<soglia<<v<<unbelt<<
            v<<"file://"<<pathbase.toStdString()<<"module_map.jpeg"<<
            v<<"file://"<<pathbase.toStdString()<<"plot_IF.jpg"<<
            v<<"file://"<<pathbase.toStdString()<<"plot_corrections.jpg"<<
            v<<"file://"<<pathbase.toStdString()<<"plot_alpha.jpg"<<
            v<<qdtst.toStdString()<<
            v<<1900+timeinfo->tm_year<<"-"<<1+timeinfo->tm_mon<<"-"<<timeinfo->tm_mday<<'\n';
    myfile.close();

    msg="waiting for new job";
    textEdit -> append(msg);

    // release the image
    //smap.release();
    //smapSlope.release();
    waitKey(100);
}


void VISfieldR::plotMap(){
    if(iProcess==0)
        return;
    //****************  creation of the module maps IF and slopeDeviations
    printf("->plotMap()");
    smap=Mat::zeros(Nirs,Njrs,CV_8UC3);
    smapSlope=Mat::zeros(Nirs,Njrs,CV_8UC3);
    msg="    Module_map image: height = " + QString::number(Nirs) +"  width = " +QString::number(Njrs);
    textEdit -> append (msg);
    double slopeDev,intfat;
    double slopeDevThreshold=DSBslopeThres -> value();
    slopeDevThreshold=slopeDevThreshold/1000.;//rad, this the dinamic range of the map
    irsStart=0;
    if(irsStep==-1){
        irsStart=Nirs;
    }
    //printf("irsStart=%d irsStep=%d\n",irsStart,irsStep);
    int iSmap=irsStart;
    for(int irs=0;irs<=Nirs;irs++){
        for(int jrs=0;jrs<Njrs;jrs++){
            if(mappa.at<double>(irs,jrs,0)>1.){
                intfat=mappa.at<double>(irs,jrs,3);
                setRGB(intfat);
                smap.at<Vec3b>(iSmap,jrs)=vRGB;
                slopeDev=mappa.at<double>(irs,jrs,2);
                slopeDev=(slopeDev+slopeDevThreshold)/2./slopeDevThreshold;
                setRGB(slopeDev);
                smapSlope.at<Vec3b>(iSmap,jrs)=vRGB;
            }
        }
        iSmap=iSmap+irsStep;
    }
    // show Intercept-factor map (smap)
    imshow("Module map", smap );
    int Nrow=smap.rows;
    int Ncol=smap.cols;
    resizeWindow("Module map",Ncol,Nrow);
    //show smapSlope
    imshow("Module Slope-Deviation map",smapSlope);
    resizeWindow("Module Slope-Deviation map",Ncol,Nrow);
    waitKey(100);

    // save smap
    QString filesave;
    filesave=pathfile+"module_map.jpeg";
    imwrite(filesave.toStdString().c_str(), smap );
    filesave=pathfile+"module_mapSlope.jpeg";
    imwrite(filesave.toStdString().c_str(), smapSlope );
    // save raw data
    QFile fileIF(pathfile+"module_map.txt");
    if (!fileIF.open(QIODevice::WriteOnly | QIODevice::Text)){
        msg="Can not open "+pathfile+"module_map.txt";
        textEdit -> append (msg);
        return;
    }
    QFile fileDS(pathfile+"module_mapSlope.txt");
    if (!fileDS.open(QIODevice::WriteOnly | QIODevice::Text)){
        msg="Can not open "+pathfile+"module_mapSlope.txt";
        textEdit -> append (msg);
        return;
    }
    QTextStream streamIF(&fileIF);
    QTextStream streamDS(&fileDS);
    iSmap=irsStart;
    for(int irs=0;irs<Nirs;irs++){
        for(int jrs=0;jrs<Njrs;jrs++){
            if(mappa.at<double>(irs,jrs,0)>1.){
                streamIF << mappa.at<double>(iSmap,jrs,3);
                streamDS << mappa.at<double>(iSmap,jrs,2);
            }
            else {
                streamIF << "NA";
                streamDS << "NA";
            }
            if(jrs<Njrs-1){
                streamIF << "\t";
                streamDS << "\t";
            }
            else {
                streamIF << "\n";
                streamDS << "\n";
            }
        }
        iSmap=iSmap+irsStep;
    }
    fileIF.close();
    fileDS.close();
    pB_crop_IF -> setEnabled(true);
    pB_crop_devSlope -> setEnabled(true);
}


void VISfieldR::crop(){
    QString msg="Please select the region to crop";
    textEdit -> append(msg);
    drag=0;
    while(drag!=-1){
        waitKey(10);
        if(drag==2) {
            Rect myROI(mpxX1,mpxY1,(mpxX2-mpxX1),(mpxY2-mpxY1));//ROI
            Cimg = imgDisplayed(myROI);// cropped image
            imshow( "ROI",Cimg);
            cropInfo[0]=1;
            cropInfo[1]=mpxX1;
            cropInfo[2]=mpxY1;
            cropInfo[3]=mpxX2-mpxX1;
            cropInfo[4]=mpxY2-mpxY1;
            waitKey(10);
            drag=-1;
        }
    }
    msg=" ";
    textEdit -> append(msg);
}


void VISfieldR::cropIF(){
    QString msg="Please select the region to crop of the Int-Fta map";
    textEdit -> append(msg);
    drag=0;
    while(drag!=-1){
        waitKey(10);
        if(drag==2) {
            Rect myROI(mpxX1,mpxY1,(mpxX2-mpxX1),(mpxY2-mpxY1));//ROI
            Cimg = smap(myROI);// cropped image
            imshow( "ROI",Cimg);
            cropInfo[0]=1;
            cropInfo[1]=mpxX1;
            cropInfo[2]=mpxY1;
            cropInfo[3]=mpxX2-mpxX1;
            cropInfo[4]=mpxY2-mpxY1;
            waitKey(10);
            drag=-1;
        }
    }
    if(mpxX1>mpxX2){
        int iTmp=mpxX1;
        mpxX1=mpxX2;
        mpxX2=iTmp;
    }
    if(mpxY1>mpxY2){
        int jTmp=mpxY1;
        mpxY1=mpxY2;
        mpxY2=jTmp;
    }
    int cont=0,iSmap;
    double sumx=0.;
    for(int jrs=mpxX1;jrs<=mpxX2;jrs++){
        for(int irs=mpxY1;irs<=mpxY2;irs++){
            iSmap=irs;
            if(irsStep==-1)
                iSmap=Nirs-irs;
            if(mappa.at<double>(iSmap,jrs,0)>0.){
                sumx=sumx+mappa.at<double>(iSmap,jrs,3);
                cont++;
            }
        }
    }
    msg="<Int-Fat>="+QString::number(sumx/cont)+" nData="+QString::number(cont);
    textEdit -> append(msg);
}


void VISfieldR::cropDevSlope(){
    QString msg="Please select the region to crop of the Int-Fta map";
    textEdit -> append(msg);
    drag=0;
    while(drag!=-1){
        waitKey(10);
        if(drag==2) {
            Rect myROI(mpxX1,mpxY1,(mpxX2-mpxX1),(mpxY2-mpxY1));//ROI
            Cimg = smapSlope(myROI);// cropped image
            imshow( "ROI",Cimg);
            cropInfo[0]=1;
            cropInfo[1]=mpxX1;
            cropInfo[2]=mpxY1;
            cropInfo[3]=mpxX2-mpxX1;
            cropInfo[4]=mpxY2-mpxY1;
            waitKey(10);
            drag=-1;
        }
    }
    if(mpxX1>mpxX2){
        int iTmp=mpxX1;
        mpxX1=mpxX2;
        mpxX2=iTmp;
    }
    if(mpxY1>mpxY2){
        int jTmp=mpxY1;
        mpxY1=mpxY2;
        mpxY2=jTmp;
    }
    int cont=0,iSmap;
    double sumx=0.;
    for(int jrs=mpxX1;jrs<=mpxX2;jrs++){
        for(int irs=mpxY1;irs<=mpxY2;irs++){
            iSmap=irs;
            if(irsStep==-1)
                iSmap=Nirs-irs;
            if(mappa.at<double>(iSmap,jrs,0)>0.){
                sumx=sumx+mappa.at<double>(iSmap,jrs,2);
                cont++;
            }
        }
    }
    msg="<devSlope>(mrad)="+QString::number(sumx/cont*1000.)+" nData="+QString::number(cont);
    textEdit -> append(msg);
}


void VISfieldR::setPoint(int Npoint){
    int iROImax,jROImax;
    Rect myROI(mpxX1,mpxY1,(mpxX2-mpxX1),(mpxY2-mpxY1));//ROI
    Mat croppedImage = imgDisplayed(myROI);// cropped image
    jROImin=mpxX1;
    iROImin=mpxY1;
    jROImax=mpxX2;
    iROImax=mpxY2;
    imshow( "ROI",croppedImage);
    waitKey(10);
    QString msg="click on the POI with medium button";
    textEdit -> append(msg);
    do{
        waitKey(100);
    } while(drag!=10);
    msg="press LeftButton to confirm";
    textEdit -> append(msg);
    do{
        mpxX1=sB_j->value();
        mpxY1=sB_i->value();
        line(croppedImage,Point(mpxX1-jROImin,0),Point(mpxX1-jROImin,iROImax-1),Scalar(0,255,0), 1);
        line(croppedImage,Point(0,mpxY1-iROImin),Point(jROImax-1,mpxY1-iROImin),Scalar(0,255,0), 1);
        imshow( "ROI",croppedImage);
        imshow("Image",imgDisplayed);
        waitKey(10);
        imgDisplayed=img.clone();
        croppedImage = imgDisplayed(myROI);
    }while(drag!=1);
    mpxX1=sB_j->value();
    mpxY1=sB_i->value();
    px[Npoint][0]=mpxX1-cx;
    px[Npoint][1]=mpxY1-cy;
    printf("->setPoint(%d): j=%d i=%d -> px[%d][0]= %f\tpx[%d][1]= %f \n",mpxX1,mpxY1,Npoint,Npoint,px[Npoint][0],Npoint,px[Npoint][1]);
    croppedImage.release();
    drag=3;
}


void VISfieldR::setJPp(){
    Np=0;
    iSetRef=0;
    viewFrameX();
    textEdit -> append(msgPoint[0]);
    drag=0;
    do{
        waitKey(10);
        if(drag==2)
            setPoint(0);
    }while(drag!=3);
    textEdit -> append(msgPoint[1]);
    drag=0;
    do{
        waitKey(10);
        if(drag==2)
            setPoint(1);
    }while(drag!=3);
    if(Np==0)
        Np=2;
    textEdit -> append(msgPoint[8]);
}



void VISfieldR::setCamera(){
    printf("setCamera() Nt=%d Np=%d...\n",Nt,Np);
    iSetRef=0;//camera has to be set
    int NtUsed=Nt;
    QMessageBox msgBox1;
    QString msg="";
    if(Np<2){
        msgBox1.setText("ATTENTION!");
        msgBox1.setInformativeText("Before you must set the joint-pivot points\n");
        msgBox1.setStandardButtons(QMessageBox::Ok);
        msgBox1.setDefaultButton(QMessageBox::Ok);
        msgBox1.exec();
        return;
    }
    Qt::CheckState state2;
    state2 =  checkBox_x0here->checkState();
    if( state2 == Qt::Checked ){
        NtUsed=Nt-2;
        // msgBox1.setText("ATTENTION!");
        // msgBox1.setInformativeText("Do you want to select another image?\n");
        // msgBox1.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        // msgBox1.setDefaultButton(QMessageBox::No);
        // int ret = msgBox1.exec();
        // if(ret==QMessageBox::Yes)
        //     imgNameVertex = QFileDialog::getOpenFileName(this,tr("Choose the Img get at x=0"),pathfile);
        // lineEdit_img ->setText(imgNameVertex.section('/', -4));
        // QString last=imgNameVertex.section('/',-1);
        Nimg=sB_Nimg -> value();
        NimgRef=Nimg;
        sB_NimgRef -> setValue(NimgRef);
        imgNameVertex=pathfile+bNf+QString::number(Nimg).rightJustified(4, '0')+eImg;
        valT=JPGtime(imgNameVertex);
        showframe(imgNameVertex);
    }
    cropInfo[0]=0;//no visualization of the stored ROI
    Pc[0]=  dSB_Xc -> value();
    Pc[1]=  dSB_Yc -> value();
    Pc[2]=  dSB_Zc -> value();
    EA[0]=  dSB_Yaw -> value(); //yaw
    EA[1]=  dSB_Pitch -> value(); //pitch
    EA[2]=  dSB_Roll -> value(); //roll
    for(int i=0;i<3;i++){
        xSto[i]=Pc[i];
        EA[i]= EA[i]/180.*pig;
        xSto[i+3]=EA[i];
    }
    Qt::CheckState state;
    state =  checkBox_fullLength->checkState();
    double v2jp=dSB_vjp ->value();
    double yp=0.25/focal*xmax*xmax;
    double iFullL=1.0;
    if( state == Qt::Checked )
        iFullL=0.0;
    // //lato destro attuato nuovo modulo ENI
    // P[0][0]=0.;
    // P[0][1]=6550.;//length/2.+940.;//7365.;//length/2.-461;//<<<<<<<<<<<461 mm = rientro target adesivi lato destro!!!!
    // P[0][2]=385.;
    // // //lato sinistro NON attuato nuovo modulo ENI croce su traliccio
    // // P[1][0]=0;
    // // P[1][1]=-6101.;//-length/2-100.;//-6490.;//-length/2.+461;//<<<<<<<<<<<461 mm = rientro target adesivi lato sinistro!!!!
    // // P[1][2]=245.;
    // //lato sinistro NON attuato nuovo modulo ENI croce su flangia
    // P[1][0]=0;
    // P[1][1]=-5890.;//-length/2-100.;//-6490.;//-length/2.+461;//<<<<<<<<<<<461 mm = rientro target adesivi lato sinistro!!!!
    // P[1][2]=0.;

    //markers for joint-pivot axis
    P[0][0]=doubleSpinBox_xR -> value();
    P[0][1]=doubleSpinBox_yR -> value();
    P[0][2]=doubleSpinBox_zR -> value();
    P[1][0]=doubleSpinBox_xL -> value();
    P[1][1]=doubleSpinBox_yL -> value();
    P[1][2]=doubleSpinBox_zL -> value();

    //corner DOWN RIGHT
    P[2][0]=xmax;
    P[2][1]=lengthDW/2.-iFullL*(panel_width+gapDW/2.);
    //corner TOP RIGHT
    P[3][0]=-xmax;
    P[3][1]=lengthUP/2.-iFullL*(panel_width+gapUP/2.);
    //corner TOP LEFT
    P[4][0]=-xmax;
    P[4][1]=-(lengthUP/2.-iFullL*(panel_width+gapUP/2.));
    //corner DOWN LEFT
    P[5][0]=+xmax;
    P[5][1]=-(lengthDW/2.-iFullL*(panel_width+gapDW/2.));
    for(int i=2;i<6;i++)
        P[i][2]=yp-v2jp;
    P[6][0]=0;
    P[6][1]=length/2.;
    P[6][2]=focal-v2jp;
    P[7][0]=0;
    P[7][1]=-length/2.;
    P[7][2]=focal-v2jp;
    for(int i=0;i<NtUsed;i++)
        printf("P%d      %f      %f        %f\n", i, P[i][0],P[i][1],P[i][2]);
    imgDisplayed=img.clone();
    imshow("Image",imgDisplayed);
    QString filetmp=pathbase+"imgDisplayed.JPG";
    imwrite(filetmp.toStdString().c_str(),imgDisplayed);
    waitKey(10);
    Qt::CheckState state1;
    state1 =  checkBox_useOldPoints->checkState();
    if( state1 == Qt::Unchecked ){
        Np=2;//Npoints
        textEdit -> append(msgPoint[Np]);
        drag=0;
        while(drag!=-1 && Np<NtUsed){
            waitKey(10);
            if(drag==2) {
                setPoint(Np);
                Np++;
                textEdit -> append(msgPoint[Np]);
            }
            if(drag==-1){
                printf("target N=%d: non visibile\n",Np);
                px[Np][0]=1.0e+9;
                px[Np][1]=1.0e+9;
                Np++;
            }
            drag=0;
            waitKey(10);
        }
        msg="Done! Npoints ="+QString::number(Np)+" Ntargets = "+QString::number(NtUsed);
        textEdit -> append(msg);
    }
    if(NtUsed!=Np ){
        msgBox1.setText("ATTENTION!");
        msgBox1.setInformativeText("Ntargets="+QString::number(NtUsed)+" != Npoints="+QString::number(Np)+"\n");
        msgBox1.setStandardButtons(QMessageBox::Ok);
        msgBox1.setDefaultButton(QMessageBox::Ok);
        msgBox1.exec();
    }
    // impostazione e lancio fit
    int n;
    if( state2 == Qt::Unchecked )
        n=6;
    else{
        n=5;
        Pc[0]=0.;//the camera is at x=0!!!
    }
    int m=NtUsed*2;//4 spigoli +2 fulcri
    //int m=6*2+2;//6 punti (x,y) + 2 distanze
    int lwa=m*n+5*n+m;
    int iwa[n];
    double x[n],fvec[m],wa[lwa];
    double tol=sqrt(dpmpar(1));
    if(n==6){
        x[0]=Pc[0];
        x[1]=Pc[1];
        x[2]=Pc[2];
        x[3]=EA[0];
        x[4]=EA[1];
        x[5]=EA[2];
    }
    else{
        Pc[0]=0.;
        x[0]=Pc[1];
        x[1]=Pc[2];
        x[2]=EA[0];
        x[3]=EA[1];
        x[4]=EA[2];
    }
    pTF[0].Nt=NtUsed;
    printf("BestFit with m=%d NtUsed=%d n=%d\n",m, NtUsed,n);
    int info=lmdif1(fcn, &pTF, m, n, x,fvec, tol, iwa, wa, lwa);
    if(n==6){
        Pc[0]=x[0];
        Pc[1]=x[1];
        Pc[2]=x[2];
        EA[0]=x[3];
        EA[1]=x[4];
        EA[2]=x[5];
    }
    else{
        Pc[1]=x[0];
        Pc[2]=x[1];
        EA[0]=x[2];
        EA[1]=x[3];
        EA[2]=x[4];
    }
    msg="Fit: info = "+QString::number(info);
    textEdit -> append(msg);
    lineEdit_chi2->setText(QString::number(chi2,'f',1));
     dSB_Xc -> setValue(Pc[0]);
     dSB_Yc -> setValue(Pc[1]);
     dSB_Zc -> setValue(Pc[2]);
     dSB_Yaw   -> setValue(EA[0]*180./pig); //yaw
     dSB_Pitch -> setValue(EA[1]*180./pig); //pitch
     dSB_Roll  -> setValue(EA[2]*180./pig); //roll
    for(int i=0;i<NtUsed;i++){
        printf("Target i = %d **************************\n",i);
        printf("P :    %f      %f        %f\n", P[i][0],P[i][1],P[i][2]);
        printf("px calcolato:    %f      %f\n",px[i][2],px[i][3]);
        printf("px sperimentale: %f      %f\n",px[i][0],px[i][1]);
    }
    QFile ftP(fileResult);//per salvataggio risultati
    int iOKopen=1;
    if(!ftP.open(QIODevice::ReadWrite | QIODevice::Text))
        iOKopen=0;
    QTextStream stream ( &ftP );
    //traccia dot sui target e sul punto calcolato
    img.copyTo(Cimg);
    for(int i=0;i<NtUsed;i++){
        if (px[i][0] < 0.9e+8){
            Pt1.x=int(px[i][0]+cx);
            Pt1.y=int(px[i][1]+cy);
            //circle(Cimg, Pt1, 10, cvScalar(65535), -1);
            circle(Cimg, Pt1, 10, Scalar(255), -1);
            Pt1.x=int(px[i][2]+cx)-10;
            Pt1.y=int(px[i][3]+cy)-10;
            Pt2.x=Pt1.x+20;
            Pt2.y=Pt1.y+20;
            //rectangle(Cimg, Pt1,Pt2, cvScalar(65535),1, 8, 0 );
            rectangle(Cimg, Pt1,Pt2, Scalar(255),1, 8, 0 );
        }
        if(iOKopen==1)
            stream<<px[i][0]+cx<<"\t"<<px[i][1]+cy<<"\t"<<px[i][2]+cx<<"\t"<<px[i][3]+cy<<"\n";
    }
    imshow( "ROI",Cimg);
    msgBox1.setText("ATTENTION!");
    msgBox1.setInformativeText("Do you want to use this image and data as reference?\n");
    msgBox1.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox1.setDefaultButton(QMessageBox::Yes);
    int ret = msgBox1.exec();
    if(ret==QMessageBox::Yes){
        iSetRef=1;
        for(int i=0;i<300;i++)
            IntFact[i][0]=0.;
        if(iOKopen==1){
            for(int i=0;i<n;i++)
                stream<<"\t"<<x[i];
            stream<<"\n";
        }
        if( state2 == Qt::Unchecked ){
            //storage Nframe position and attitude
            NimgRef=Nimg;
            if(iOKopen==1){
                stream << Nimg<<"\n";
                ftP.close();
            }
            sB_NimgRef ->setValue(Nimg);
            sB_Nimg->setValue(Nimg);
            viewFrameX();
        }else{
            if(iOKopen==1){
                stream << imgNameVertex<<"\n";
                ftP.close();
            }
            sB_NimgRef ->setValue(NimgRef);
            thetaImg=0.;
            showframe(imgNameVertex);
        }
        sB_Nimg -> setValue(NimgRef);
    }
    else{
        iSetRef=0;
         dSB_Xc -> setValue(xSto[0]);
         dSB_Yc -> setValue(xSto[1]);
         dSB_Zc -> setValue(xSto[2]);
         dSB_Yaw   -> setValue(xSto[0]*180./pig); //yaw
         dSB_Pitch -> setValue(xSto[1]*180./pig); //pitch
         dSB_Roll  -> setValue(xSto[2]*180./pig); //roll
         viewFrameX();
    }
    if(iOKopen==1) ftP.close();
}


void VISfieldR::pointer(int i1B2G3R, double thetaI, double x1, double y1, double z1){
    double v[3];
    Mod2Dro(0.,thetaI,0.,x1,y1,z1);
    for(int i=0;i<3;i++)
        v[i]=Vservice[i];
    xyz2pxpy(v[0],v[1],v[2]);//compute pxX pxY given X,Y,Z,Pc[3],EA[3]
    Pt1.x=int(pxX+cx);
    Pt1.y=int(pxY+cy);
    if(i1B2G3R==1)
        circle(imgDisplayed, Pt1, 10, Scalar(255,0,0), -1);
    else if(i1B2G3R==2)
        circle(imgDisplayed, Pt1, 10, Scalar(0,255,0), -1);
    else if(i1B2G3R==3)
        circle(imgDisplayed, Pt1, 10, Scalar(0,0,255), -1);
}


void VISfieldR::liner(int i1B2G3R, double thetaI, double x1, double y1, double z1,
                      double x2, double y2, double z2){
    double v[3];
    Mod2Dro(0.,thetaI,0.,x1,y1,z1);
    for(int i=0;i<3;i++)
        v[i]=Vservice[i];
    xyz2pxpy(v[0],v[1],v[2]);
    Pt1.x=int(pxX+cx);//j
    Pt1.y=int(pxY+cy);//i
    Mod2Dro(0.,thetaI,0.,x2,y2,z2);
    for(int i=0;i<3;i++)
        v[i]=Vservice[i];
    xyz2pxpy(v[0],v[1],v[2]);
    Pt2.x=int(pxX+cx);
    Pt2.y=int(pxY+cy);
    if(i1B2G3R==1)
        line(imgDisplayed,Pt1,Pt2, Scalar(255,0,0), 3);
    else if(i1B2G3R==2)
        line(imgDisplayed,Pt1,Pt2, Scalar(0,255,0), 3);
    else if(i1B2G3R==3)
        line(imgDisplayed,Pt1,Pt2, Scalar(0,0,255), 3);
}



void VISfieldR::ROI(double thetaI,double xTop, double xBot ){
    double x1,y1,z1,v[3],Lpx1,Lpx2;
    double Plim[4][2];
    double v2jp=dSB_vjp ->value();
    printf("->ROI(thetaI=%f,xTop=%f,xBot=%f)\n",thetaI,xTop,xBot);
//left Top
    x1=xTop;
    if(x1<0)
        y1=-lengthUP/2.;
    else
        y1=-lengthDW/2.;
    z1=0.25/focal*x1*x1-v2jp;
    Mod2Dro(0.,thetaI,0.,x1,y1,z1);
    for(int i=0;i<3;i++)
        v[i]=Vservice[i];
    xyz2pxpy(v[0],v[1],v[2]);
    Plim[0][1]=pxX+cx;//j
    Plim[0][0]=pxY+cy;//i
    //printf("leftTopji: %f %f\n",Plim[0][1],Plim[0][0]);
//right Top
    y1=-y1;
    Mod2Dro(0.,thetaI,0.,x1,y1,z1);
    for(int i=0;i<3;i++)
        v[i]=Vservice[i];
    xyz2pxpy(v[0],v[1],v[2]);
    Plim[1][1]=pxX+cx;
    Plim[1][0]=pxY+cy;
    //printf("rightTopji: %f %f\n",Plim[1][1],Plim[1][0]);
    Lpx1=sqrt((Plim[1][1]-Plim[0][1])*(Plim[1][1]-Plim[0][1])+
            (Plim[1][0]-Plim[0][0])*(Plim[1][0]-Plim[0][0]));
    dpximg=y1/Lpx1;
//left Bot
    x1=xBot;
    if(x1<0)
        y1=-lengthUP/2.;
    else
        y1=-lengthDW/2.;
    z1=0.25/focal*x1*x1-v2jp;
    Mod2Dro(0.,thetaI,0.,x1,y1,z1);
    for(int i=0;i<3;i++)
        v[i]=Vservice[i];
    xyz2pxpy(v[0],v[1],v[2]);
    Plim[2][1]=pxX+cx;
    Plim[2][0]=pxY+cy;
    //printf("leftBotji: %f %f\n",Plim[2][1],Plim[2][0]);
//right Bot
    y1=-y1;
    Mod2Dro(0.,thetaI,0.,x1,y1,z1);
    for(int i=0;i<3;i++)
        v[i]=Vservice[i];
    xyz2pxpy(v[0],v[1],v[2]);
    Plim[3][1]=pxX+cx;
    Plim[3][0]=pxY+cy;
    //printf("rightBotji: %f %f\n",Plim[3][1],Plim[3][0]);
    Lpx2=sqrt((Plim[2][1]-Plim[3][1])*(Plim[2][1]-Plim[3][1])+
            (Plim[2][0]-Plim[3][0])*(Plim[2][0]-Plim[3][0]));
    dpximg=dpximg+y1/Lpx2;
    Lpx=(Lpx1+Lpx2)/2.;
    //calcolo coefficienti angolari e termine noto
    if(fabs(Plim[1][1]-Plim[0][1]) > 0.000001)
        mjt=(Plim[1][0]-Plim[0][0])/(Plim[1][1]-Plim[0][1]);//limite superiore
    else
        mjt=0.;
    cjt=-mjt*Plim[0][1]+Plim[0][0];
    if(fabs(Plim[3][1]-Plim[2][1]) > 0.000001)
        mjb=(Plim[3][0]-Plim[2][0])/(Plim[3][1]-Plim[2][1]);//limite inferiore
    else
        mjb=0.;
    cjb=-mjb*Plim[2][1]+Plim[2][0];
    if(fabs(Plim[2][0]-Plim[0][0]) > 0.000001)
        mil=(Plim[2][1]-Plim[0][1])/(Plim[2][0]-Plim[0][0]);//limite sinistro
    else
        mil=0.;
    cil=-mil*Plim[0][0]+Plim[0][1];
    if(fabs(Plim[3][0]-Plim[1][0]) > 0.000001)
        mir=(Plim[3][1]-Plim[1][1])/(Plim[3][0]-Plim[1][0]);//limite destro
    else
        mir=0.;
    cir=-mir*Plim[1][0]+Plim[1][1];
    jLeft=(Plim[0][1]+Plim[2][1])/2.;
    iLeft=(Plim[0][0]+Plim[2][0])/2.;
    jRight=(Plim[1][1]+Plim[3][1])/2.;
    iRight=(Plim[1][0]+Plim[3][0])/2.;
    // printf("ROI(%f, %f, %f):\nmjt=%f cjt=%f\n mjb=%f cjb=%f \nmil=%f cil=%f\n mir=%f cir=%f\ndpximg=%f\n",
    //        thetaI,xTop,xBot,  mjt,   cjt,     mjb,   cjb,     mil,   cil,     mir,   cir, dpximg);
    // printf("\tjLeft=%f iLeft=%f Lpx=%f\n",jLeft,iLeft,Lpx);
}



void VISfieldR::GnuPlot(string script){//Create plot by GnuPlot and script command
    string command="gnuplot -persist "+pathbase.toStdString()+script;
    printf("GnuPlot with command: %s\n",command.c_str());
    int ierr=system(command.c_str());
    if(ierr!=0) {
        QString msg="GnuPlot returns error for command= "+QString::fromStdString(command);
        QMessageBox msgBox;
        msgBox.setText(msg);
        msgBox.exec();
    }
}




int fcn(void *p, int m, int n, const double *x, double *fvec, int iflag){
    /* calculate the functions at x and return the values in fvec[0] through fvec[m-1] */
    struct pointToFit *pTF = (struct pointToFit *)p;
    int fredeg=m-n;
    int Ntt=pTF[0].Nt;
    int k=0;
    chi2=0.;
    if(n==6){
        Pc[0]=x[0];
        Pc[1]=x[1];
        Pc[2]=x[2];
        EA[0]=x[3];
        EA[1]=x[4];
        EA[2]=x[5];
    }
    else{
        Pc[1]=x[0];
        Pc[2]=x[1];
        EA[0]=x[2];
        EA[1]=x[3];
        EA[2]=x[4];
    }
    for(int i=0;i<Ntt;i++){
        if(px[i][0] < 0.9e+9){
            xyz2pxpy(P[i][0],P[i][1],P[i][2]);//compute pxX pxY given X,Y,Z,Pc[3],EA[3]
//                if(i == 0){
//                 printf("\n>>>>>>>>>>>>>>>>>> drone data:\n");
//                 printf("Pc:    %f      %f        %f\n",x[0],x[1],x[2]);
//                 printf("EA:    %f      %f        %f\n",x[3]*180./pig,x[4]*180./pig,x[5]*180./pig);
//                }
//            printf("Target i = %d **************************\n",i);
//            printf("P :    %f      %f        %f\n", P[i][0],P[i][1],P[i][2]);
//            printf("px calcolato:    %f      %f\n",pxX,pxY);
//            printf("px sperimentale: %f      %f\n",px[i][0],px[i][1]);
            px[i][2]=pxX;
            px[i][3]=pxY;
            if(i<6){
                fvec[k]=pxX-px[i][0];
                chi2=chi2+fvec[k]*fvec[k];
                //printf("        fvec@%d&%d: %f        ",k+1,k+2,fvec[k]);
                k++;
                fvec[k]=pxY-px[i][1];
                chi2=chi2+fvec[k]*fvec[k];
                //printf(" %f\n",fvec[k]);
                //printf("chi2= %f\n",chi2);
                k++;
            }
            else{
                i++;
                xyz2pxpy(P[i][0],P[i][1],P[i][2]);
                px[i][2]=pxX;
                px[i][3]=pxY;
//                printf("Target i = %d **************************\n",i);
//                printf("P :    %f      %f        %f\n", P[i][0],P[i][1],P[i][2]);
//                printf("px calcolato:    %f      %f\n",pxX,pxY);
//                printf("px sperimentale: %f      %f\n",px[i][0],px[i][1]);
                double mr=(pxY-px[i-1][3])/(pxX-px[i-1][2]);// y=mr*x+qr
                double qr=px[i-1][3]-mr*px[i-1][2];
//                printf("mr=%f\tqr=%f\n",mr,qr);
                fvec[k]=(px[i-1][1]-(mr*px[i-1][0]+qr))/sqrt(1.+mr*mr);
                chi2=chi2+fvec[k]*fvec[k];
                k++;
                fvec[k]=(px[i][1]-(mr*px[i][0]+qr))/sqrt(1.+mr*mr);
                chi2=chi2+fvec[k]*fvec[k];
//                printf("chi2= %f\n",chi2);
            }
        }
    }
     // printf("k=%d\n",k);
     //printf(">>> chi2 = %f\n\tPc: %f ; %f ; %f\n\tEA: %f ; %f ; %f\n",
     //       chi2,x[0],x[1],x[2],x[3],x[4],x[5]);
    chi2fin=chi2/double(fredeg);
    return(0);
}



void unitV12(double x1, double y1, double z1, double x2, double y2, double z2){
    // calculates the unit vector P1->P2
    double sum=0.;
    Vservice[0]=x2-x1;
    Vservice[1]=y2-y1;
    Vservice[2]=z2-z1;
    for(int i=0;i<3;i++)
        sum=sum+Vservice[i]*Vservice[i];
    sum=sqrt(sum);
    for(int i=0;i<3;i++)
        Vservice[i]=Vservice[i]/sum;
}



void xyz2pxpy(double X, double Y, double Z){ //compute pxX pxY given X,Y,Z,Pc[3],EA[3]
    double v[3],vr[3];
    unitV12(Pc[0],Pc[1],Pc[2],X,Y,Z);//versore drone -> target
    v[0]=Vservice[0];
    v[1]=Vservice[1];
    v[2]=Vservice[2];
    Mod2Dro(EA[0],EA[1],EA[2],v[0],v[1],v[2]);//trasformazione di v nel rif. del drone
    vr[0]=Vservice[0];
    vr[1]=Vservice[1];
    vr[2]=Vservice[2];
    pxX=vr[0]*fx/vr[2];
    pxY=vr[1]*fy/vr[2];
}


void Mod2Dro(double yaw,double pitch, double roll,double xm, double ym, double zm){
    double c1,s1,c2,s2,c3,s3;
    c1=cos(roll);
    s1=sin(roll);
    c2=cos(pitch);
    s2=sin(pitch);
    c3=cos(yaw);
    s3=sin(yaw);
    // Vservice[0]=             c2*c3*xm               +c2*s3*ym           -s2*zm;
    // Vservice[1]=( s1*s2*c3-c1*s3 )*xm    +(s1*s2*s3+c1*c3)*ym        +s1*c2*zm;
    // Vservice[2]=( c1*s2*c3+s1*s3 )*xm    +(c1*s2*s3-s1*c3)*ym        +c1*c2*zm;
    //Earth2Plane (VSIproPT)
    Vservice[0]=           c2*c3*xm            -c2*s3*ym      +s2*zm;
    Vservice[1]=(c1*s3+c3*s1*s2)*xm +(c1*c3-s1*s2*s3)*ym   -c2*s1*zm;
    Vservice[2]=(s1*s3-c1*c3*s2)*xm +(c3*s1+c1*s2*s3)*ym   +c1*c2*zm;
}



void imgDimension(double tg2p){
  int sign,signOld;
  double xx=fabs(xP);
  double xr=xx,step,ze,zeOld,tetaMIN,tetaMAX,xzero;
  double fOr=(focal+0.25/focal*xr*xr-df)/radius;//<<<<<<<ERRORE f0r=focal/radius
  double ctphi;
  if(tg2p*(1.0-fOr*fOr)+1.0 >= 0.0)
   ctphi= (fOr*tg2p+sqrt(tg2p*(1.0-fOr*fOr)+1.0))/(1.0+tg2p);
  else
   ctphi=1.0/fOr;
  double alpha0 = 2.0*atan(0.5/focal*xx);
  tetaMIN = min(pig-alpha0+acos(ctphi),pig-alpha0-acos(ctphi));
  tetaMAX = max(pig-alpha0+acos(ctphi),pig-alpha0-acos(ctphi));
  printf("->imgDimension: tg2p=%f fOr=%f xx=%f tphi=%f alpha0=%f tetaMIN=%f tetaMAX=%f\n",tg2p,fOr,xx,ctphi,alpha0,tetaMIN,tetaMAX);

  for(int ic=1;ic<=4;ic++){
    if(ic < 3)
      step=5.0;
    else
      step=-5.0;
    if(ic == 1 || ic==3) xr=xx;
    do{
     zeOld=zero(ic,tetaMIN,tetaMAX,xr);
     xr=xr+step;
     ze=zero(ic,tetaMIN,tetaMAX,xr);
     if(abs(ze)>abs(zeOld)) step=-step;
//     printf("xr = %f\n",xr);
//     waitKey(0);
    } while(abs(ze)>abs(zeOld));
    sign=segno(ze);
    signOld=segno(zeOld);
    while(sign*signOld > 0.0){
      signOld=sign;
      zeOld=ze;
      xr=xr+step;
      ze = zero(ic,tetaMIN,tetaMAX,xr);
      sign=segno(ze);
//      printf("    ic=%d xr=%f ze=%f\n",ic,xr,ze);
//      waitKey(0);
    }
    xzero=(xr-step)-zeOld*step/(ze-zeOld);
    //printf("\tic=%d step=%f xzero=%f\n",ic,step,xzero);
    if(ic==1){
      if(xP>0.)
        xsMax=xzero;
      else
        xsMin=-xzero;
    }
    else if(ic==2){
      if(xP>0.)
        xrMax=xzero;
      else
        xrMin=-xzero;
    }
    else if(ic==3){
      if(xP>0.)
        xsMin=xzero;
      else
        xsMax=-xzero;
    }
    else if(ic==4){
      if(xP>0.)
        xrMin=xzero;
      else
        xrMax=-xzero;
    }
  }
//  double xBari=0.5*xrMax+0.5*xrMin;
  if(xsMin<xrMin) xsMin=xrMin;
  if(xsMax>xrMax) xsMax=xrMax;
  printf("\txP=%f dview=%f xrMin=%f xsMin=%f xsMax=%f xrMax=%f\n",xP,dview,xrMin,xsMin,xsMax,xrMax);
  Magne=(xrMax-xrMin)/2.0/radius;
//  MagneTheory= (dview-y)/(y+focal);
//  printf("    Magnitude_exp = %f  Magnitude_theory = %f ratio = %f\n",Magne,MagneTheory,Magne/MagneTheory);
}


double zero(int ic, double tetaMIN, double tetaMAX, double xr){
      double yr,ct,ai,an,ar,ze,xp=0.0,yp=0.0;
      yr = 0.25*(xr*xr)/focal;
      if(ic == 1){
        xp=radius*sin(tetaMIN);//thetaMAX
        yp=radius*cos(tetaMIN)+focal-df;//<<<<<<<--------modifica
      }
      else if(ic==2){
        ct = (-radius*(focal-df-yr)+
          sqrt(-4.0*(focal-df)*yr*radius*radius+xr*xr*(focal-df+yr)*(focal-df+yr)))/(focal-df+yr)/(focal-df+yr);//<<<<<<<--------modifica
        xp=radius*sqrt(1.0-ct*ct);
        yp=radius*ct+focal-df;//<<<<<<<--------modifica
      }
      else if(ic==3){
        xp=radius*sin(tetaMAX);
        yp=radius*cos(tetaMAX)+focal-df;//<<<<<<<--------modifica
      }
      else if(ic==4){
        ct = (-radius*(focal-df-yr)-
          sqrt(-4.0*(focal-df)*yr*radius*radius+xr*xr*(focal-df+yr)*(focal-df+yr)))/(focal-df+yr)/(focal-df+yr);//<<<<<<<--------modifica
        xp=-radius*sqrt(1.0-ct*ct);
        yp=radius*ct+focal-df;//<<<<<<<--------modifica
      }
      ai = atan((yp-yr)/(xr-xp));
      if(ai<0.0) ai=ai+pig;
      an = atan(2.0*focal/xr);
      if(an<0.0) an=an+pig;
      ar = pig/2.0+atan((fabs(xP)-xr)/(dview-yr));
      ze = an-0.5*(ai+ar);
      return(ze);
}



int segno(double val){
 if (val > 0.) return 1;
 if (val < 0.) return -1;
return 0;
}



void setRGB(double x){
    uchar Blue=0;
    uchar Green=0;
    uchar Red=0;
    if(x>=0.0 && x< 0.25){
        Blue=255;
        Green=uchar(255.*4.*x+0.5);
        Red=0;
    }
    else if(x>= 0.25 && x< 0.50){
        Blue=uchar(255.*(1.-4.*(x-0.25)));
        Green=255;
        Red=0;
    }
    else if(x>= 0.50 && x< 0.75){
        Blue=0;
        Green=255;
        Red=uchar(255.*4.*(x-0.50));
    }
    else if(x>= 0.75 && x< 1.){
        Blue=0;
        Green=uchar(255.*(1.-4.*(x-0.75)));
        Red=255;
    }
    else if(x>= 1.){
        Blue=255;
        Green=255;
        Red=255;
    }
    vRGB.val[0] = Blue;
    vRGB.val[1] = Green;
    vRGB.val[2] = Red;
    //cout<<Blue<<"\t"<<Green<<"\t"<<Red<<"\n";
}
