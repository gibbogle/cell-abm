// myvtk.cpp

#include <vtkCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>

#ifdef _WIN32
#include "windows.h"
#endif
#include "myvtk.h"
#include "log.h"
#include "transfer.h"

LOG_USE();

// Define interaction style
class MouseInteractorStyle4 : public vtkInteractorStyleTrackballCamera
{
  public:
	static MouseInteractorStyle4* New();
	vtkTypeMacro(MouseInteractorStyle4, vtkInteractorStyleTrackballCamera);

	virtual void OnLeftButtonDown()
	{
	  leftb = true;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
	}

	virtual void OnMiddleButtonDown()
	{
//	  std::cout << "Pressed middle mouse button." << std::endl;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
	}

	virtual void OnRightButtonDown()
	{
//	  std::cout << "Pressed right mouse button." << std::endl;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnRightButtonDown();
	}

	virtual void OnLeftButtonUp()
	{
//	  std::cout << "Released left mouse button." << std::endl;
//	  LOG_QMSG("Released left mouse button.");
	  leftb = false;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
	}

	virtual void OnMiddleButtonUp()
	{
//	  std::cout << "Released middle mouse button." << std::endl;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnMiddleButtonUp();
	}

	virtual void OnRightButtonUp()
	{
//	  std::cout << "Released right mouse button." << std::endl;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnRightButtonUp();
	}

};

vtkStandardNewMacro(MouseInteractorStyle4);

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
MyVTK::MyVTK(QWidget *page, QWidget *key_page)
{
	zoomlevel = 0.7;
	double backgroundColor[] = {0.0,0.0,0.0};


	Pi = 4*atan(1.0);
    page_VTK = page;
	leftb = false;
    key_canvas(key_page);
    qvtkWidget = new QVTKWidget(page,QFlag(0));
	LOG_MSG("Created a new QVTKWidget");
	QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(qvtkWidget);

	// Associate the layout with page_VTK
    page_VTK->setLayout(layout);

	// Create a renderer, and add it to qvtkWidget's render window.
	// The renderer renders into the render window. 
	ren = vtkRenderer::New();     
    renWin = qvtkWidget->GetRenderWindow();
    renWin->AddRenderer(ren);
	ren->SetBackground(backgroundColor);
	ren->ResetCamera();
	iren = qvtkWidget->GetInteractor();

	vtkSmartPointer<MouseInteractorStyle4> style = vtkSmartPointer<MouseInteractorStyle4>::New();
	iren->SetInteractorStyle( style );

	iren->Initialize();

	// Create mappers
	createMappers();

	// Create image filter for save Snapshot()
//	w2img = vtkWindowToImageFilter::New();
//	pngwriter = vtkSmartPointer<vtkPNGWriter>::New();
//	jpgwriter = vtkSmartPointer<vtkJPEGWriter>::New();

	first_VTK = true;
	DCmotion = false;
    DCfade = false;
	playing = false;
	paused = false;
    opacity = 1.0;
    display_celltype[1] = true;
    display_celltype[2] = true;
    TCpos_list.clear();
    ren->GetActiveCamera()->Zoom(zoomlevel);		// try zooming OUT
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
MyVTK::~MyVTK()
{
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::key_canvas(QWidget *key_page)
{
//    QGraphicsScene* scene = new QGraphicsScene(QRect(0, 0, 130, 280));
    QGraphicsScene* scene = new QGraphicsScene(QRect(0, 0, 130, 310));
    QBrush brush;
    QGraphicsTextItem *text;

//    brush.setColor(QColor(255,128,77));
    brush.setColor(QColor(150,100,0));      // dark brown
    brush.setStyle(Qt::SolidPattern);
	scene->addEllipse(10,10,20,20,Qt::NoPen, brush);
	text = scene->addText("FDC");
    text->setPos(35, 10);

//    brush.setColor(QColor(255,77,128));
    brush.setColor(QColor(200,60,100));     // dark red
    brush.setStyle(Qt::SolidPattern);
    scene->addEllipse(10,40,20,20,Qt::NoPen, brush);
    text = scene->addText("MRC");
    text->setPos(35, 40);

    brush.setColor(QColor(30,20,255));      // dark blue
    scene->addEllipse(10,70,20,20,Qt::NoPen, brush);
	text = scene->addText("Naive B cell");
    text->setPos(35, 70);

    brush.setColor(QColor(0,200,255));      // light blue
    scene->addEllipse(10,100,20,20,Qt::NoPen, brush);
	text = scene->addText("CCR7 UP");
    text->setPos(35, 100);

    brush.setColor(QColor(50,255,150));     // light green
    scene->addEllipse(10,130,20,20,Qt::NoPen, brush);
	text = scene->addText("EBI2 UP");
    text->setPos(35, 130);

//    brush.setColor(QColor(255,255,0));      // yellow
    brush.setColor(Qt::yellow );      // yellow
    scene->addEllipse(10,160,20,20,Qt::NoPen, brush);
	text = scene->addText("BCL6 HI");
    text->setPos(35, 160);

    brush.setColor(QColor(0,150,0));        // dark green
    scene->addEllipse(10,190,20,20,Qt::NoPen, brush);
	text = scene->addText("BCL6 LO");
    text->setPos(35, 190);

    brush.setColor(QColor(128,128,128));    // grey
    scene->addEllipse(10,220,20,20,Qt::NoPen, brush);
	text = scene->addText("Max divisions");
    text->setPos(35, 220);

    brush.setColor(QColor(255,0,0));        // red
    scene->addEllipse(10,250,20,20,Qt::NoPen, brush);
	text = scene->addText("Plasma cell");
    text->setPos(35, 250);

    brush.setColor(QColor(255,130,0));      // orange
    scene->addEllipse(10,280,20,20,Qt::NoPen, brush);
    text = scene->addText("CD4 T cell");
    text->setPos(35, 280);

	QGraphicsView* view = new QGraphicsView(key_page);
    view->setScene(scene);
//    view->setGeometry(QRect(0, 0, 150, 300));
    view->setGeometry(QRect(0, 0, 150, 330));
    view->show();
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::createMappers()
{
    vtkSphereSource *Tcell = vtkSphereSource::New();
    Tcell->SetThetaResolution(12);
    Tcell->SetPhiResolution(12);
    Tcell->SetRadius(0.5);
    TcellMapper = vtkPolyDataMapper::New();

    TcellMapper->SetInputConnection(Tcell->GetOutputPort());

//	vtkSphereSource *Dcell = vtkSphereSource::New();
//	Dcell->SetThetaResolution(12);
//	Dcell->SetPhiResolution(12);
//	Dcell->SetRadius(1.0);
//	DcellMapper = vtkPolyDataMapper::New();
//	DcellMapper->SetInputConnection(Dcell->GetOutputPort());
//	vtkCylinderSource *bond = vtkCylinderSource::New();
//	bond->SetResolution(12);
//	bond->SetRadius(0.15);
//	bond->SetHeight(1);
//	bondMapper = vtkPolyDataMapper::New();
//	bondMapper->SetInputConnection(bond->GetOutputPort());

	double rSphere0 = 0.5;
	double rSphere1 = 0.3;
	double rCylinder = 0.2;
	// create sphere geometry
	vtkSphereSource *sphere0 = vtkSphereSource::New();
	sphere0->SetRadius(rSphere0);
	sphere0->SetThetaResolution(18);
	sphere0->SetPhiResolution(18);
	sphere0->SetCenter(0.0,0.0,0.0);
	vtkPolyData *sData0 = sphere0->GetOutput();

	vtkSphereSource *sphere1 = vtkSphereSource::New();
	sphere1->SetRadius(rSphere1);
	sphere1->SetThetaResolution(18);
	sphere1->SetPhiResolution(18);
	sphere1->SetCenter(0.0,-1.0,0.0);
	vtkPolyData *sData1 = sphere1->GetOutput();

	vtkSphereSource *sphere2 = vtkSphereSource::New();
	sphere2->SetRadius(rSphere1);
	sphere2->SetThetaResolution(18);
	sphere2->SetPhiResolution(18);
	sphere2->SetCenter(0.0,1.0,0.0);
	vtkPolyData *sData2 = sphere2->GetOutput();

	// create cylinder geometry
	vtkCylinderSource *cylinder = vtkCylinderSource::New();
	cylinder->SetCenter(0.0, 0.0, 0.0);
	cylinder->SetRadius(rCylinder);
	cylinder->SetHeight(2.0);
	cylinder->SetResolution(18);
	vtkPolyData *cData = cylinder->GetOutput();

	// Append the data
	vtkAppendPolyData* append1 = vtkAppendPolyData::New();
	append1->AddInput(cData);
	append1->AddInput(sData1);
	append1->AddInput(sData2);

	vtkPolyData *dumbell1 = append1->GetOutput();

	vtkTransform *t2 = vtkTransform::New();
	t2->PostMultiply();
	t2->RotateZ(90);
	vtkTransformPolyDataFilter *tf2 = vtkTransformPolyDataFilter::New();
	tf2->SetTransform(t2);
	tf2->SetInput(dumbell1);
	vtkPolyData *dumbell2 = tf2->GetOutput();

	vtkTransform *t3 = vtkTransform::New();
	t3->PostMultiply();
	t3->RotateX(90);
	vtkTransformPolyDataFilter *tf3 = vtkTransformPolyDataFilter::New();
	tf3->SetTransform(t3);
	tf3->SetInput(dumbell1);
	vtkPolyData *dumbell3 = tf3->GetOutput();

	vtkAppendPolyData* append2 = vtkAppendPolyData::New();
	append2->AddInput(sData0);
	append2->AddInput(dumbell1);
	append2->AddInput(dumbell2);
	append2->AddInput(dumbell3);

	// Rendering objects.
//	FDcellMapper = vtkPolyDataMapper::New();
//	FDcellMapper->SetInput(append2->GetOutput());

	// Is this OK?
	sphere0->Delete();
	sphere1->Delete();
	sphere2->Delete();
	append1->Delete();
	append2->Delete();
}


//-----------------------------------------------------------------------------------------
// The cell info is fetched from the DLL by ExecThread::snapshot().
// The info is transmitted in the integer arrays cell_list[] and DC_list[]
// The info is transferred here into TCpos_list and DCpos_list, which are Qlists.
//-----------------------------------------------------------------------------------------
void MyVTK::get_cell_positions(bool fast)
{
    int ninfo = 6;
//    LOG_QMSG("get_cell_positions");
	double BC_diam = 0.9;
//	double DC_diam = 1.8;
    TCpos_list.clear();
//	DCpos_list.clear();
//	bondpos_list.clear();
    for (int i=0; i<ncell_list; i++) {
        int j = ninfo*i;
		CELL_POS cp;
        cp.tag = cell_list[j];
        cp.x = cell_list[j+1];
        cp.y = cell_list[j+2];
        cp.z = cell_list[j+3];
        cp.state = cell_list[j+4];
        cp.diameter = cell_list[j+5]/100.;
//        cp.diameter = BC_diam;
        TCpos_list.append(cp);
//        double r, g, b;
//        if (cp.state > 0) {
//            unpack(cp.state, &r, &g, &b);
//        } else {
//            r = g = b = 0;
//        }
//        sprintf(msg,"B cell: %d: tag: %d pos: %d %d %d state: %d %lf %lf %lf",i,cp.tag,cp.x,cp.y,cp.z,cp.state,r,g,b);
//        LOG_MSG(msg);
	}

}

//-----------------------------------------------------------------------------------------
// Reading cell positions from a file.  Out of date.
// NOT USED
//-----------------------------------------------------------------------------------------
void MyVTK::read_cell_positions(QString infileName, QString outfileName, bool savepos)
{	
    TCpos_list.clear();
//    DCpos_list.clear();
//    bondpos_list.clear();
	QString line, saveline;
	QTextStream *out = NULL;
	QFile *vtkdata = NULL;
	if (savepos) {
		vtkdata = new QFile(outfileName);
		if (!vtkdata->open(QFile::Append )) {
			LOG_MSG("Open failure on vtk file");
			return;
		}
		out = new QTextStream(vtkdata);
	}
	QFile posdata(infileName);
	if (posdata.open(QFile::ReadOnly)) {
		QTextStream in(&posdata);
		do {
			line = in.readLine();
			if (line.length() > 0) {
				if (savepos) {
					*out << line << "\n";
					out->flush();
				}
				saveline = line;
				QStringList s = line.split(" ",QString::SkipEmptyParts);
				if (s[0].compare("T") == 0) {
					CELL_POS cp;
					cp.tag = s[1].toInt();
					cp.x = s[2].toInt();
					cp.y = s[3].toInt();
					cp.z = s[4].toInt();
					cp.diameter = s[5].toDouble();
					cp.state = s[6].toDouble();
                    TCpos_list.append(cp);
//				}
//                else if (s[0].compare("D") == 0) {
//					CELL_POS cp;
//					cp.tag = s[1].toInt();
//					cp.x = s[2].toInt();
//					cp.y = s[3].toInt();
//					cp.z = s[4].toInt();
//					cp.diameter = s[5].toDouble();
//					cp.state = s[6].toDouble();
//					DCpos_list.append(cp);
//				} else if (s[0].compare("B") == 0) {
//					BOND_POS cp;
//					cp.BCtag = s[1].toInt();
//					cp.DCtag = s[2].toInt();
//					bondpos_list.append(cp);
				} else if (s[0].compare("E") == 0) {
					break;
				} 
			}
		} while (!line.isNull());
	}
	posdata.close();
	if (savepos) {
		delete out;
		vtkdata->close();
		delete vtkdata;
	}

	if (QFile::exists(infileName)) {
		QFile::rename(infileName,"TO_REMOVE");
	}
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::init()
{
    T_Actor_list.clear();
//	D_Actor_list.clear();
//	Bnd_Actor_list.clear();
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::cleanup()
{
	int i;
//	vtkActor *actor;
    ACTOR_TYPE a;

	LOG_MSG("VTK cleanup");
    for (i = 0; i<T_Actor_list.length(); i++) {
        a = T_Actor_list[i];
        ren->RemoveActor(a.actor);
	}
//	for (i = 0; i<D_Actor_list.length(); i++) {
//        a = D_Actor_list[i];
//        ren->RemoveActor(a.actor);
//	}
//	for (i = 0; i<Bnd_Actor_list.length(); i++) {
//		actor = Bnd_Actor_list[i];
//        ren->RemoveActor(actor);
//	}
    T_Actor_list.clear();
//	D_Actor_list.clear();
//	Bnd_Actor_list.clear();
	first_VTK = true;	
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::renderCells(bool redo, bool zzz)
{
    process_Tcells();
	if (first_VTK) {
		LOG_MSG("Initializing the renderer");
		ren->ResetCamera();
	}
	iren->Render();
    first_VTK = false;
//    if (record) {
//        recorder();
//    }
}

//---------------------------------------------------------------------------------------------
// Interprets an int as rgb
// USE_CELLTYPE_COLOUR, the cell type is passed in cp.state, and the colours are those
// that were chosen in the GUI.
// Otherwise cp.state is packed (r,g,b)
//---------------------------------------------------------------------------------------------
void MyVTK::unpack(int x, double *rr, double *gg, double *bb)
{
	int z, r, g, b;

    if (USE_CELLTYPE_COLOUR) {

    } else {
        z = x;
        r = z>>16;
        z = r;
        z = z<<16;

        x = x - z;

        z = x;
        g = z>>8;
        z = g;
        z = z<<8;

        b = x - z;
    }
    *rr = r/255.;
    *gg = g/255.;
    *bb = b/255.;
}

//---------------------------------------------------------------------------------------------
// This improved, simpler method assumes only that the cell tag (cp.tag) is unique.
// This tag is the index into T_Actor_list[], which grows as the maximum tag increases.
// T_Actor_list now includes the field .active, which is maintained in sync with the
// status of .actor in the renderer object ren, i.e. when the actor is added to the scene
// .active is set to true, and when the actor is removed .active is set to false.
// In each scene update step the status of a tag is recorded in in_pos_list[].
// in_pos_list is compared with T_Actor_list, and actors corresponding to tags that have been
// dropped (i.e. have in_pos_list[tag] = false but A_Actor_list[tag].active = true) are marked
// as .active = false in T_Actor_list, and removed from the scene.
//---------------------------------------------------------------------------------------------
void MyVTK::process_Tcells()
{
    int i, tag, maxtag;
    double r, g, b, h;
    double growth;
	CELL_POS cp;
	int axis_centre = -2;	// identifies the ellipsoid centre
	int axis_end    = -3;	// identifies the ellipsoid extent in 5 directions
	int axis_bottom = -4;	// identifies the ellipsoid extent in the -Y direction, i.e. bottom surface
    bool dbug = false;
    ACTOR_TYPE a;
    ACTOR_TYPE *ap;
    QColor qcolor, colour1, colour2;
    bool display;

 //   LOG_QMSG("process_Tcells");
    int np = TCpos_list.length();
    if (np == 0) return;
    int na = T_Actor_list.length();
    if (istep < 0) {
        sprintf(msg,"na: %d np: %d",na,np);
        LOG_MSG(msg);
        dbug = true;
    }
    colour1 = celltype_colour[1];
    colour2 = celltype_colour[2];
    maxtag = 0;
	for (i=0; i<np; i++) {
        cp = TCpos_list[i];
        tag = cp.tag;
        maxtag = max(tag,maxtag);
	}
    if (dbug) {
        sprintf(msg,"maxtag: %d",maxtag);
        LOG_MSG(msg);
    }
    ap = &a;
    for (tag=na; tag<=maxtag; tag++) {
        ap->actor = vtkActor::New();
        ap->actor->SetMapper(TcellMapper);
        ap->active = false;
        T_Actor_list.append(a);
    }
    na = T_Actor_list.length();
    bool *in_pos_list;
    in_pos_list = new bool[na];
    for (tag=0; tag<na; tag++)
        in_pos_list[tag] = false;
    for (i=0; i<np; i++) {
        cp = TCpos_list[i];
        tag = cp.tag;
        growth = cp.state/100.;
        display = true;
        // Wnt_Cyclin
        int itype;
        if (growth == 0) {          // G0
            r = 0.5;
            g = 0.5;
            b = 0;
            itype = 0;
        } else if (growth < 0.4) {         // G1
            r = 1;
            g = 1;
            b = 0;
            itype = 1;
        } else if (growth < 0.6) {  // S
            r = 1;
            g = 0.5;
            b = 0;
            itype = 2;
        } else if (growth < 0.9) {  // G2
            r = 1;
            g = 0.5;
            b = 0;
            itype = 3;
        } else {                    // M
            r = 1;
            g = 0.5;
            b = 0;
            itype = 4;
        }
        if (USE_CELLTYPE_COLOUR) {
            qcolor = celltype_colour[itype];
            r = qcolor.red()/255.;
            g = qcolor.green()/255.;
            b = qcolor.blue()/255.;
            display = display_celltype[itype];
        }
        if (!display) continue;

        in_pos_list[tag] = true;
        if (dbug) {
            sprintf(msg,"i: %d tag: %d",i,tag);
            LOG_MSG(msg);
        }
        ap = &T_Actor_list[tag];
        if (!ap->active) {  // Make active an actor with new tag in TCpos_list
            if (dbug) {
                sprintf(msg,"adding actor: %d",tag);
                LOG_MSG(msg);
            }
            ren->AddActor(ap->actor);
            ap->active = true;
        }
//        if (USE_CELLTYPE_COLOUR) {
//            qcolor = celltype_colour[cp.state];
//            r = qcolor.red()/255.;
//            g = qcolor.green()/255.;
//            b = qcolor.blue()/255.;
////                sprintf(msg,"celltype: %d  r,g,b: %f %f %f",cp.state,r,g,b);
////                LOG_MSG(msg);
//        } else {
//            unpack(cp.state, &r, &g, &b);
//        }


        /*
        // Cell_cycle
        if (growth < 1.51) {        // G1
            r = 0;
            g = 1;
            b = 0;
            qcolor = colour1;
        } else if (growth < 1.92) { // S
            r = 1;
            g = 0.5;
            b = 0;
            qcolor = colour2;
        } else if (growth < 2.11) { // G2
            r = 1;
            g = 0.5;
            b = 0;
            qcolor = colour2;
        } else {                    // M
            r = 1;
            g = 0.5;
            b = 0;
            qcolor = colour2;
        }


        // Wnt_Cyclin
        if (growth < 0.4) {         // G1
            r = 1;
            g = 1;
            b = 0;
        } else {                    // S, G2, M
            r = 1;
            g = 0.5;
            b = 0;
//        } else if (growth < 0.6) {  // S
//            r = 1;
//            g = 0.5;
//            b = 0;
//        } else if (growth < 0.9) {  // G2
//            r = 1;
//            g = 0.5;
//            b = 0;
//        } else {                    // M
//            r = 1;
//            g = 0.5;
//            b = 0;
        }
        */

        ap->actor->GetProperty()->SetColor(r, g, b);
        ap->actor->SetScale(cp.diameter);
//        if (cp.highlight == 0)
//            ap->actor->GetProperty()->SetOpacity(opacity);
        ap->actor->SetPosition(cp.x, cp.y, cp.z);
	}
    for (int k=0; k<T_Actor_list.length(); k++) {
        ap = &T_Actor_list[k];
        if (ap->active && !in_pos_list[k]) {
            if (dbug) {
                sprintf(msg,"removing actor: %d",k);
                LOG_MSG(msg);
            }
            ren->RemoveActor(ap->actor);
            ap->active = false;
        }
	}
}


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool MyVTK::startPlayer(QString posfile, QTimer *theTimer, bool save)
{
	save_image = save;
	LOG_QMSG(posfile);
	timer = theTimer;
	playerData = new QFile(posfile);
	if (!playerData->open(QFile::ReadOnly)) {
		LOG_MSG("Open failure on VTK file");
		return false;
	}
	playerStream = new QTextStream(playerData);
	if (!first_VTK) {
		cleanup();
	}
	playing = true;
	paused = false;

	if (save_image) {
        w2i = vtkWindowToImageFilter::New();
        w2i->SetInput(renWin);	//the render window
//		writer = vtkSmartPointer<vtkPNGWriter>::New();
        jpgwriter = vtkSmartPointer<vtkJPEGWriter>::New();
        jpgwriter->SetInputConnection(w2i->GetOutputPort());
		framenum = 0;
		LOG_MSG("set up writer");
	}
	LOG_MSG("playing");
	return true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool MyVTK::nextFrame()
{
	LOG_MSG("VTK: nextFrame");
	if (!playing)
		return false;
	if (paused)
		return true;
	if (playerStream->atEnd()) {
		LOG_MSG("nextFrame: no more data");
		stop();
		return false;
	}
    TCpos_list.clear();
//	DCpos_list.clear();
//	bondpos_list.clear();
	int k = 0;
	QString line;
	do {
		line = playerStream->readLine();
		if (line.length() > 0) {
			k++;
			QStringList s = line.split(" ",QString::SkipEmptyParts);
			if (s[0].compare("T") == 0) {
				CELL_POS cp;
				cp.tag = s[1].toInt();
				cp.x = s[2].toInt();
				cp.y = s[3].toInt();
				cp.z = s[4].toInt();
				cp.diameter = s[5].toDouble();
				cp.state = s[6].toDouble();
                TCpos_list.append(cp);
//			} else if (s[0].compare("D") == 0) {
//				CELL_POS cp;
//				cp.tag = s[1].toInt();
//				cp.x = s[2].toInt();
//				cp.y = s[3].toInt();
//				cp.z = s[4].toInt();
//				cp.diameter = s[5].toDouble();
//				cp.state = s[6].toDouble();
//				DCpos_list.append(cp);
//			} else if (s[0].compare("B") == 0) {
//				BOND_POS cp;
//				cp.BCtag = s[1].toInt();
//				cp.DCtag = s[2].toInt();
//				bondpos_list.append(cp);
			} else if (s[0].compare("E") == 0) {
				break;
			}
		}
	} while (true);

	bool redo = false;
	if (first_VTK) {
		redo = true;
	}
    renderCells(redo,false);
	char numstr[5];
	sprintf(numstr,"%04d",framenum);
	if (save_image) {
        w2i->Modified();	//important
        jpgwriter->SetFileName((casename + numstr + ".jpg").toStdString().c_str());
        jpgwriter->Write();
	}
	framenum++;
	return true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::pause()
{
    paused = true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::playon()
{
    paused = false;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::stop()
{
    if (save_image) {
        if (jpgwriter) jpgwriter->Delete();
        if (pngwriter) pngwriter->Delete();
        w2i->Delete();
    }
    delete playerStream;
    playerData->close();
    delete playerData;
    timer->stop();
    playing = false;
    paused = false;
}


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::saveSnapshot(QString fileName, QString imgType)
{
    vtkSmartPointer<vtkWindowToImageFilter> w2img = vtkWindowToImageFilter::New();

	w2img->SetInput(renWin);
	if (imgType.compare("png") == 0) {
		vtkSmartPointer<vtkPNGWriter> pngwriter = vtkPNGWriter::New();
		pngwriter->SetInputConnection(w2img->GetOutputPort()); 
		w2img->Modified();
		pngwriter->SetFileName((fileName.toStdString()).c_str()); 
		pngwriter->Write();
//		pngwriter->Delete();	// Note: using vtkSmartPointer, delete is not necessary.
	} else if (imgType.compare("jpg") == 0) {
		vtkJPEGWriter *jpgwriter = vtkJPEGWriter::New();
		jpgwriter->SetInputConnection(w2img->GetOutputPort()); 
		w2img->Modified();
		jpgwriter->SetFileName((fileName.toStdString()).c_str()); 
		jpgwriter->Write();
//		jpgwriter->Delete();
	} else if (imgType.compare("tif") == 0) {
		vtkTIFFWriter *tifwriter = vtkTIFFWriter::New();
		tifwriter->SetInputConnection(w2img->GetOutputPort()); 
		w2img->Modified();
		tifwriter->SetFileName((fileName.toStdString()).c_str()); 
		tifwriter->Write();
//		tifwriter->Delete();
	} else if (imgType.compare("bmp") == 0) {
		vtkBMPWriter *bmpwriter = vtkBMPWriter::New();
		bmpwriter->SetInputConnection(w2img->GetOutputPort()); 
		w2img->Modified();
		bmpwriter->SetFileName((fileName.toStdString()).c_str()); 
		bmpwriter->Write();
	}
}


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::set_celltype_colour(COLOUR_TYPE *colour, QString str)
{
    if (str.compare("red") == 0) {
        colour->r = 1.0;
        colour->g = 0.0;
        colour->b = 0.0;
    } else if (str.compare("orange") == 0) {
        colour->r = 0.8;
        colour->g = 0.5;
        colour->b = 0.0;
    } else if (str.compare("yellow") == 0) {
        colour->r = 1.0;
        colour->g = 1.0;
        colour->b = 0.0;
    } else if (str.compare("green") == 0) {
        colour->r = 0.0;
        colour->g = 1.0;
        colour->b = 0.0;
    } else if (str.compare("blue") == 0) {
        colour->r = 0.0;
        colour->g = 0.0;
        colour->b = 1.0;
    } else if (str.compare("purple") == 0) {
        colour->r = 1.0;
        colour->g = 0.0;
        colour->b = 1.0;
    } else if (str.compare("brown") == 0) {
        colour->r = 0.5;
        colour->g = 0.5;
        colour->b = 0.2;
    }

}

/*
//-----------------------------------------------------------------------------------------
// Uses: w2i, renWin, pngwriter,
//       record, record_basename, record_nframes, record_it, framenum
//-----------------------------------------------------------------------------------------
void MyVTK::startRecorder(QString basefile, int nframes)
{
    if (w2i == 0) {
        w2i = vtkWindowToImageFilter::New();
        w2i->SetInput(renWin);	//the render window
//        jpgwriter = vtkSmartPointer<vtkJPEGWriter>::New();
//        jpgwriter->SetInputConnection(w2i->GetOutputPort());
        pngwriter = vtkSmartPointer<vtkPNGWriter>::New();
        pngwriter->SetInputConnection(w2i->GetOutputPort());
}
    framenum = 0;
    LOG_MSG("set up writer");
    record = true;
    record_basename = basefile;
    record_nframes = nframes;
    record_it = 0;

    LOG_MSG("Started recording");
}

//-----------------------------------------------------------------------------------------
// Uses: w2i, videoOutput, pngwriter,
//       record, record_basename, record_nframes, record_it, framenum
//-----------------------------------------------------------------------------------------
void MyVTK::recorder()
{
    char numstr[5];
    char filename[512];

    sprintf(msg,"recorder: record_it: %d",record_it);
    LOG_MSG(msg);
    if (record_it > record_nframes) {
        record = false;
        stopRecorder();
        return;
    }
    vtkImageData *id = vtkImageData::New();
    id = w2i->GetOutput();
    w2i->Modified();	//important
    id->Update();
    int width = id->GetDimensions()[0];
    int height = id->GetDimensions()[1];
    if (width == 0) {
        LOG_QMSG("ERROR: recorder: vtkImageData dimension = 0");
        exit(1);
    }
    if (RECORD_VIDEO) {
        if (!videoOutput->isOpen()) {
            // Generate temporary filename
            tempFile = new QTemporaryFile("qt_temp.XXXXXX.avi");
            if (tempFile->open())
            {
               // Open media file and prepare for recording
               QString fileName = tempFile->fileName();
                bool recording = videoOutput->openMediaFile(width, height, fileName.toAscii().data());
                if (!recording) {
                    LOG_QMSG("ERROR: openMediaFile failed");
                    record = false;
                    return;
                }
            }
        }
        bool success = videoOutput->newVtkFrame(id);
        if (!success) {
            LOG_QMSG("ERROR: newVtkFrame failed");
            record = false;
            exit(1);
        }
        record_it++;
        return;
    }
//    strcpy(filename,record_basename);
//    strcat(filename,numstr);
//    strcat(filename,".jpg");
//    strcpy(filename ,(record_basename + numstr + ".jpg").toStdString().c_str());
//    jpgwriter->SetFileName(filename);
//    jpgwriter->Write();
    sprintf(numstr,"%05d",framenum);
    strcpy(filename ,(record_basename + numstr + ".png").toStdString().c_str());
    pngwriter->SetFileName(filename);
    pngwriter->Write();
    sprintf(msg,"recorder: it: %d frame: %d filename: %s  id dimensions: %d %d",record_it,framenum,filename,id->GetDimensions()[0],id->GetDimensions()[1]);
    LOG_MSG(msg);
    framenum++;
    record_it++;
}

//-----------------------------------------------------------------------------------------
// Uses: record
//-----------------------------------------------------------------------------------------
void MyVTK::stopRecorder()
{
    record = false;
    if (RECORD_VIDEO) {
        videoOutput->closeMediaFile();
//        QString fileName = QFileDialog::getSaveFileName(page_VTK,
//                                                        "Save File",    //tr("Save File"),
//                                                        QString(),
//                                                        "Videos (*.avi)");  //tr("Videos (*.avi)"));
        QString fileName = "testfile.avi";
        if (fileName.isNull() == false) {
           QFile::copy(tempFile->fileName(), fileName);
        }
        delete tempFile;
        tempFile = 0x0;
    }
//    jpgwriter->RemoveAllInputs();
//    jpgwriter->Delete();
//    w2i->RemoveAllInputs();
//    w2i->Delete();
    LOG_MSG("Stopped recording");
}
*/


/*
void MyVTK::recordingSlot()
{
   if (recording)
   {
      videoOutput->closeMediaFile();
      recordingButton->setIcon(QIcon(":/robotnavigator/resources/video-camera-png.png"));
      recordingButton->setToolTip("Record Video");
      recordingTimer->stop();
      recording = false;
      QString fileName = QFileDialog::getSaveFileName(this,
                                                      tr("Save File"),
                                                      QString(),
                                                      tr("Videos (*.avi)"));
      if (fileName.isNull() == false)
      {
         QFile::copy(tempFile->fileName(), fileName);
      }
      delete tempFile;
      tempFile = 0x0;
   }
   else
   {
      blinkCount = 0;
      // Generate temporary filename
      tempFile = new QTemporaryFile("qt_temp.XXXXXX.avi");
      if (tempFile->open())
      {
         // Open media file and prepare for recording
         QString fileName = tempFile->fileName();
         recording = videoOutput->openMediaFile(640, 480, fileName.toAscii().data());
      }
      // Change tool tip on video button and start blinking timer
      recordingButton->setToolTip("Stop Recording");
      recordingTimer->start(500);
   }
}
*/

/*
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::process_Dcells()
{
    int i, tag, maxtag;
    double r, g, b;
    CELL_POS cp;
    vtkActor *actor;
//    double FDCColor[] = {1.0,0.5,0.3};
    double FDCColor[] = {0.59,0.39,0.0};  // 150,100,0
//    double MRCColor[] = {1.0,0.3,0.5};
    double MRCColor[] = {0.78,0.24,0.39};     // 200,60,100
    bool dbug = false;
    ACTOR_TYPE a;
    ACTOR_TYPE *ap;

//    LOG_QMSG("process_Dcells");
    int np = DCpos_list.length();
    int na = D_Actor_list.length();
    if (istep < 0) {
        sprintf(msg,"na: %d np: %d",na,np);
        LOG_MSG(msg);
        dbug = true;
    }
    maxtag = 0;
    for (i=0; i<np; i++) {
        cp = DCpos_list[i];
        tag = cp.tag;
        maxtag = max(tag,maxtag);
    }
    if (dbug) {
        sprintf(msg,"maxtag: %d",maxtag);
        LOG_MSG(msg);
    }
    ap = &a;
    for (tag=na; tag<=maxtag; tag++) {
        ap->actor = vtkActor::New();
        ap->actor->SetMapper(FDcellMapper);
        ap->active = false;
        D_Actor_list.append(a);
    }
    na = D_Actor_list.length();
    bool *in_pos_list;
    in_pos_list = new bool[na];
    for (tag=0; tag<na; tag++)
        in_pos_list[tag] = false;
    for (i=0; i<np; i++) {
        cp = DCpos_list[i];
        tag = cp.tag;
        in_pos_list[tag] = true;
        if (dbug) {
            sprintf(msg,"i: %d tag: %d",i,tag);
            LOG_MSG(msg);
        }
        ap = &D_Actor_list[tag];
        if (!ap->active) {  // Make active an actor with new tag in TCpos_list
            if (dbug) {
                sprintf(msg,"adding actor: %d",tag);
                LOG_MSG(msg);
            }
            ren->AddActor(ap->actor);
            ap->active = true;
        }
        if (cp.state == 100)
            ap->actor->GetProperty()->SetColor(FDCColor);
        else if (cp.state == 200)
            ap->actor->GetProperty()->SetColor(MRCColor);
        ap->actor->SetPosition(cp.x, cp.y, cp.z);
    }
    for (int k=0; k<D_Actor_list.length(); k++) {
        ap = &D_Actor_list[k];
        if (ap->active && !in_pos_list[k]) {
            if (dbug) {
                sprintf(msg,"removing actor: %d",k);
                LOG_MSG(msg);
            }
            ren->RemoveActor(ap->actor);
            ap->active = false;
        }
    }
}


//---------------------------------------------------------------------------------------------
// A cylinder is created orientated along the y-axis, i.e. along b = (0,1,0)
// To change the orientation to the vector v, we first create a vector r
// normal to both b and v: r = bxv, this will be the axis of rotation.
// We now need to rotate the cylinder by theta about r, where theta is the angle between
// b and v, i.e. sin(theta) = |r|/(|b||v|) = |r|/|v|
// We can now use actor.RotateWXYZ(theta,r[0],r[1],r[2]) where theta is in degrees
// What is bxv when b = (0,1,0) and v = (v0,v1,v2)?
// r = [v[2],0,-v[0]]
//---------------------------------------------------------------------------------------------
void MyVTK::process_bonds()
{
    int i, j;
    BOND_POS bp;
    vtkActor *actor, *B_actor, *D_actor;
    double bpos[3], v[3];
    double Pi = 3.15159;
    double *bcpos, *dcpos;
    double bondColor[] = {0.5,0.0,0.0};

    int na = Bnd_Actor_list.length();
    int np = bondpos_list.length();

    // First remove all old bonds (strictly speaking we should remove only those not in the new list)

    for (int k=0; k<na; k++) {
        ren->RemoveActor(Bnd_Actor_list[k]);
    }

    Bnd_Actor_list.clear();

    for (i=0; i<np; i++) {
        bp = bondpos_list[i];
        actor = vtkActor::New();
        actor->SetMapper(bondMapper);
        actor->GetProperty()->SetColor(bondColor);
        B_actor = Bnd_Actor_list[bp.BCtag];
        if (B_actor != 0)
            bcpos = B_actor->GetPosition();
        else {
            sprintf(msg,"B_actor = 0 in bond: %d %d",i,bp.BCtag);
            LOG_MSG(msg);
            exit(1);
        }
        D_actor = D_Actor_list[bp.DCtag];
        if (D_actor != 0)
            dcpos = D_actor->GetPosition();
        else {
            sprintf(msg,"D_actor = 0 in bond: %d %d",i,bp.DCtag);
            LOG_MSG(msg);
            exit(1);
        }

        for (j=0; j<3; j++) {
            bpos[j] = (bcpos[j] + dcpos[j])/2;
            v[j] = bcpos[j] - dcpos[j];
        }
        double v_mod = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        double s[] = {1, v_mod, 1};
        actor->SetScale(s);
        for (j=0; j<3; j++)
            v[j] = v[j]/v_mod;

        double sina = sqrt(v[0]*v[0] + v[2]*v[2]);
        double cosa = v[1];
        double theta = asin(sina)*(180.0/Pi);
        if (cosa < 0)
            theta = 180 - theta;

        actor->SetPosition(bpos);
        actor->RotateWXYZ(theta,v[2],0,-v[0]);
        ren->AddActor(actor);
        Bnd_Actor_list.append(actor);
    }
}
*/

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//void MyVTK::set_display_celltype(int celltype, bool display)
//{
//    display_celltype[celltype] = display;
//}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//void MyVTK::set_celltype_colour(int celltype, const QString &text)
//{
//    celltype_colour[celltype] = text;
//}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//void MyVTK::on_checkBox_CELLDISPLAY_1_toggled(bool display)
//{
//    display_celltype[1] = display;
//    renderCells(false,false);
//    LOG_QMSG("toggled display_celltype[1]");
//}
