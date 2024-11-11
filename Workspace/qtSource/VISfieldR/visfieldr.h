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

#ifndef VISFIELDR_H
#define VISFIELDR_H

#include "ui_visfieldr.h"

using namespace std;

class VISfieldR : public QWidget, private Ui::VISfieldR_DLG
{
    Q_OBJECT

public:
    VISfieldR(QWidget *parent = nullptr);
    ~VISfieldR();
    void setWin(const std::string& _winname);

private:
    Ui::VISfieldR_DLG *ui;
    void on_mouse_internal(int ev, int x, int y);
    std::string winname;
    friend void on_mouse(int ev, int x, int y, int, void* obj);

public slots:
    void crop();
    void cropIF();
    void cropDevSlope();
    void setJPp();
    void setPoint(int Npoint);
    void setCamera();
    void setPanel();
    void setGap();
    void SetFile();
    void SetMeasFile();
    void viewFrameX();
    void showframe(QString spath);
    void pointer(int i1B2G3R, double thetaI, double x1, double y1, double z1);
    void liner(int i1B2G3R, double thetaI, double x1, double y1, double z1,
               double x2, double y2, double z2);
    void saveMeasFile();
    void ROI(double thetaI, double xTop, double xBot );
    void process();
    void GnuPlot(string script);
    void updateN();
    void viewNmin();
    void viewNmax();
    void viewImgRef();
    double JPGtime(QString imgName);
    void setFontDia();
    void enableDisable();
    void plotMap();
};

#endif
