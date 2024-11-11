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
#include <QApplication>
#include "visfieldr.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    VISfieldR *dialog = new VISfieldR;
    dialog->show();
    return app.exec();
}
