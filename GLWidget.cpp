
#include "GLWidget.h"

#include <iostream>
#include <random>
#include <math.h>

#include <QGLFormat>
#include <QSvgGenerator>

#include "SystemParams.h"
#include "VertexData.h"

#include "LineCloud.h"
#include "NanoFLANN.h"

GLWidget::GLWidget(QGLFormat format, QWidget *parent) :
    QGLWidget(format, parent),
    _isMouseDown(false),
    _zoomFactor(10.0),
    _shaderProgram(0),
    _img_width(10),
    _img_height(10),
    _iterStatus(-1),
    _imageTexture(0)
    //_drawSelPoint(false)
{
}


GLWidget::~GLWidget()
{
    if(_shaderProgram) delete _shaderProgram;
    if(_imageTexture) delete _imageTexture;
}

// http://stackoverflow.com/questions/16782746/what-is-faster-than-stdpow
// http://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
inline double fastPow(double a, double b) {
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
}

MyPoint rotate_point(float cx,float cy,float angle, MyPoint p)
{
  float s = sin(angle);
  float c = cos(angle);

  // translate point back to origin:
  p.x -= cx;
  p.y -= cy;

  // rotate point
  float xnew = p.x * c - p.y * s;
  float ynew = p.x * s + p.y * c;

  // translate point back:
  p.x = xnew + cx;
  p.y = ynew + cy;
  return p;
}

void GLWidget::initializeGL()
{
    QGLFormat glFormat = QGLWidget::format();
    if (!glFormat.sampleBuffers()) { std::cerr << "Could not enable sample buffers." << std::endl; return; }

    glShadeModel(GL_SMOOTH);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClearColor( 1.0, 1.0, 1.0, 1.0 );
    glEnable(GL_DEPTH_TEST);

    _shaderProgram = new QOpenGLShaderProgram();
    if (!_shaderProgram->addShaderFromSourceFile(QOpenGLShader::Vertex, "../OrganicLabyrinth/shader.vert"))
        { std::cerr << "Cannot load vertex shader." << std::endl; return; }

    if (!_shaderProgram->addShaderFromSourceFile(QOpenGLShader::Fragment, "../OrganicLabyrinth/shader.frag"))
        { std::cerr << "Cannot load fragment shader." << std::endl; return; }

    if ( !_shaderProgram->link() )
        { std::cerr << "Cannot link shaders." << std::endl; return; }

    _shaderProgram->bind();
    _mvpMatrixLocation = _shaderProgram->uniformLocation("mvpMatrix");
    _colorLocation = _shaderProgram->attributeLocation("vertexColor");

    /*
    /*
    std::random_device rd;
    std::mt19937 e2(rd());
    std::normal_distribution<> dist(0, SystemParams::dist_std_dev);
    */

    _e2 = std::mt19937(_rd());
    _dist = std::normal_distribution<>(0, SystemParams::dist_std_dev);


    //_img_size.setWidth( 0 );
    //_img_size.setHeight( 0 );
}


bool GLWidget::event( QEvent * event )
{
    return QGLWidget::event(event);
}


// This is an override function from Qt but I can't find its purpose
void GLWidget::resizeGL(int width, int height)
{
}

void GLWidget::paintGL()
{
    if(_iterStatus == 0)
    {
        _currentIter = 0;
        //_drawSelPoint = false;
        InitCurve();
        std::cout << "init curve done\n";
    }

    if(_iterStatus == 1)
    {
        EvolveCurve();
    }

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glViewport(0, 0, this->width(),  this->height());

    int current_width = width();
    int current_height = height();

    // Set orthographic Matrix
    QMatrix4x4 orthoMatrix;

    orthoMatrix.ortho(0.0 +  _scrollOffset.x(),
                      (float)current_width +  _scrollOffset.x(),
                      (float)current_height + _scrollOffset.y(),
                      0.0 + _scrollOffset.y(),
                      -100, 100);

    // Translate the view to the middle
    QMatrix4x4 transformMatrix;
    transformMatrix.setToIdentity();
    transformMatrix.scale(_zoomFactor);

    _shaderProgram->setUniformValue(_mvpMatrixLocation, orthoMatrix * transformMatrix);

    PaintCurve();


    if(_magnitudeVao.isCreated() && SystemParams::show_mag)
    {
        int use_color_location = _shaderProgram->uniformLocation("use_color");
        _shaderProgram->setUniformValue(use_color_location, (GLfloat)0.0);

        _magnitudeTexture->bind();
        _magnitudeVao.bind();
        glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
        _magnitudeVao.release();
        _magnitudeTexture->release();
    }

    if(_imageVao.isCreated() && SystemParams::show_image)
    {
        int use_color_location = _shaderProgram->uniformLocation("use_color");
        _shaderProgram->setUniformValue(use_color_location, (GLfloat)0.0);

        _imageTexture->bind();
        _imageVao.bind();
        glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
        _imageVao.release();
        _imageTexture->release();
    }


}

// Mouse is pressed
void GLWidget::mousePressEvent(int x, int y)
{
    _isMouseDown = true;

    double dx = x + _scrollOffset.x();
    dx /= _zoomFactor;

    double dy = y + _scrollOffset.y();
    dy /= _zoomFactor;
}

// Mouse is moved
void GLWidget::mouseMoveEvent(int x, int y)
{
    double dx = x + _scrollOffset.x();
    dx /= _zoomFactor;

    double dy = y + _scrollOffset.y();
    dy /= _zoomFactor;

    // your stuff

    this->repaint();
}


// Mouse is released
void GLWidget::mouseReleaseEvent(int x, int y)
{
    _isMouseDown = false;
    double dx = x + _scrollOffset.x();
    dx /= _zoomFactor;

    double dy = y + _scrollOffset.y();
    dy /= _zoomFactor;

    // your stuff

    this->repaint();
}

void GLWidget::mouseDoubleClick(int x, int y)
{
    double dx = x + _scrollOffset.x();
    dx /= _zoomFactor;

    double dy = y + _scrollOffset.y();
    dy /= _zoomFactor;
}

void GLWidget::HorizontalScroll(int val) { _scrollOffset.setX(val); }
void GLWidget::VerticalScroll(int val) { _scrollOffset.setY(val); }
void GLWidget::ZoomIn() { this->_zoomFactor += 0.5f; }
void GLWidget::ZoomOut() { this->_zoomFactor -= 0.5f; if(this->_zoomFactor < 0.1f) _zoomFactor = 0.1f; }


void GLWidget::InitCurve()
{
    _points.clear();

    MyPoint centerPt(this->_img_width / 2, this->_img_height / 2);

    float addValue = (M_PI * 2.0 / (float)SystemParams::circle_init_slices);
    for(float a = 0.0; a < M_PI * 2.0; a += addValue)
    {
        float xPt = centerPt.x + SystemParams::circle_radius * sin(a);
        float yPt = centerPt.y + SystemParams::circle_radius * cos(a);
        _points.push_back(MyPoint(xPt, yPt));
    }

    // resampling
    ResampleCurve();

    _iterStatus = 1;

}

void GLWidget::CreateCurveVAO()
{
    // POINTS VAO
    QVector3D vecCol = QVector3D(1.0, 0.0, 0.0);
    PreparePointsVAO(_points, &_pointsVbo, &_pointsVao, vecCol);

    // LINES VAO
    vecCol = QVector3D(0.0, 0.5, 1.0);
    std::vector<MyLine> lines;
    for(size_t a = 0; a < _points.size(); a++)
    {
        if(a < _points.size() - 1) { lines.push_back(MyLine(_points[a], _points[a+1])); }
        else { lines.push_back(MyLine(_points[a], _points[0])); }
    }
    PrepareLinesVAO(lines, &_linesVbo, &_linesVao, vecCol);
}


void GLWidget::ResampleCurve()
{
    // a really simple resampling (the version from the paper)
    std::vector<MyPoint> tempPoints;

    // hack
    _points.push_back(_points[0]);

    //float threshold2 = SystemParams::D * SystemParams::k_min;

    for(size_t a = 0; a < _points.size() - 1;)
    {
        MyPoint pt1 = _points[a];
        MyPoint pt2 = _points[a + 1];

        //pt1.age++;
        tempPoints.push_back(pt1);

        float delta_factor = 0.5f * (GetDelta(pt1.x, pt1.y) + GetDelta(pt2.x, pt2.y));

        float threshold1 = SystemParams::D * SystemParams::d_split * delta_factor;

        if(pt1.Distance(pt2) > threshold1)   // split
        {
            MyPoint newPt = pt1 + (pt2 - pt1) * 0.5;
            tempPoints.push_back(newPt);
            a += 1;
        }
        else
        {
            a += 2;
        }
    }

    // rebuild curve
    _points = std::vector<MyPoint>(tempPoints);

    // hack
    float delta_factor = 0.5f * (GetDelta(_points[_points.size() - 1].x, _points[_points.size() - 1].y) + GetDelta(_points[0].x, _points[0].y));
    float threshold1 = SystemParams::D * SystemParams::d_split * delta_factor;
    if(_points[_points.size() - 1].Distance(_points[0]) < threshold1 && _points.size() > 1)
        { _points.erase(_points.begin()); }

    // uniform resampling (more precise)
    /*
    std::vector<MyPoint> tempPoints;
    UniformResample(_points, tempPoints, SystemParams::D);
    _points.clear();
    for(size_t a = 0; a < tempPoints.size(); a++)
        { _points.push_back(tempPoints[a]); }
    */
}

MyPoint GLWidget::GetDirVector(float x, float y)
{
    MyPoint retVec(0, 0);
    if(_imageVao.isCreated())
    {
        int intx = (int)x;
        int inty = (int)y;
        if(intx > 0 && intx < _img_width && inty > 0 && inty < _img_height)
        {
            retVec = MyPoint(_gradientX[intx][inty], _gradientY[intx][inty]);
        }
    }
    return retVec;
}

float GLWidget::GetDelta(float x, float y)
{
    //return SystemParams::delta_const;
    float retVal = SystemParams::delta_const;


    if(_imageVao.isCreated())
    {
        int intx = (int)x;
        int inty = (int)y;
        float spanVal = 1.0f - SystemParams::delta_bound;

        //std::cout << spanVal << "\n";

        if(spanVal > 1e-8 && intx > 0 && intx < _img_width && inty > 0 && inty < _img_height)
        {
            float val = qGray(_imgOriginal.pixel(intx, inty));
            retVal = val / 255.0;
            retVal = (retVal * spanVal) + SystemParams::delta_bound;
        }
    }

    return retVal;

}


void GLWidget::EvolveCurve()
{
    if(_currentIter == SystemParams::max_iter - 1)
    {
        _iterStatus = -1;
    }

    //clone
    std::vector<MyPoint> tempPoints(_points);

    // 1 --- BROWNIAN MOTION
    // fix me, delta is static

    for(size_t a = 0; a < _points.size(); a++)
    {
        float l1 = GetRandomNumber();
        float l2 = GetRandomNumber();
        float randX = l1 * cos(l2);
        float randY = l1 * sin(l2);

        float delta_const_D = GetDelta(_points[a].x, _points[a].y) * SystemParams::D;
        float delta_const_D_fb = delta_const_D * SystemParams::f_b;

        float xB = randX * delta_const_D_fb;
        float yB = randY * delta_const_D_fb;

        tempPoints[a] += MyPoint(xB, yB);
    }

    // 2 --- FAIRING
    // fix me, delta is static

    for(size_t a = 0; a < _points.size(); a++)
    {
        //if(_points[a].age >= maxAge) { continue;}

        MyPoint curPt = _points[a];
        MyPoint prevPt;
        MyPoint nextPt;

        if(a == 0) { prevPt = _points[_points.size() - 1]; }
        else { prevPt = _points[a - 1]; }
        if(a < _points.size() - 1) { nextPt = _points[a + 1]; }
        else { nextPt = _points[0]; }

        float prevDelta = GetDelta(prevPt.x, prevPt.y);
        float nextDelta = GetDelta(nextPt.x, nextPt.y);
        float denominator = prevDelta + nextDelta;

        float xFactor = (prevPt.x * nextDelta + nextPt.x * prevDelta) / denominator;
        float yFactor = (prevPt.y * nextDelta + nextPt.y * prevDelta) / denominator;

        xFactor -= curPt.x;
        yFactor -= curPt.y;

        tempPoints[a] += (MyPoint(xFactor, yFactor) * SystemParams::f_f);
    }

    // 3 --- ATTRACTION - REPULSION

    GetClosestPoints();

    for(size_t a = 0; a < _points.size(); a++)
    {
        MyPoint pt = GetAttractionRepulsion(a);
        tempPoints[a] += pt;

    }

    // update
    _points = std::vector<MyPoint>(tempPoints);

    // HV (work best)
    /*
    for(size_t a = 0; a < _points.size(); a++)
    {
        MyPoint pt1(this->_img_width / 2, this->_img_height / 2);
        MyPoint pt2 = tempPoints[a];
        MyPoint dirPt = pt2 - pt1;
        //dirPt = dirPt.Norm();
        MyPoint triPt(cos(dirPt.x), sin(dirPt.y));
        triPt *= 0.015;

        _points[a] = tempPoints[a] + triPt;
    }
    */


    // anisotropy
    /*for(size_t a = 0; a < _points.size(); a++)
    {
        MyPoint magNorm = GetDirVector(_points[a].x, _points[a].y).Norm();
        _points[a] = tempPoints[a] + magNorm * 0.05;
    }
    */

    // spiral !!!
    /*
    for(size_t a = 0; a < _points.size(); a++)
    {
        // _points[a] = tempPoints[a];

        MyPoint pt1 = tempPoints[a] - MyPoint(this->_img_width / 2, this->_img_height);
        MyPoint pt2 = rotate_point(0, 0, M_PI / 3.0, pt1) * 2.5f;
        MyPoint dirPt = pt2 - pt1;
        dirPt = dirPt.Norm();
        //MyPoint triPt(cos(dirPt.x), sin(dirPt.y));
        dirPt *= 0.1;

        _points[a] = tempPoints[a] + dirPt;

    }
    */


    // BLOCKS
    /*
    float frequency = 5;
    for(size_t a = 0; a < _points.size(); a++)
    {
        // _points[a] = tempPoints[a];
        MyPoint pt = tempPoints[a];

        double fractpartx, fractparty, intpartx, intparty;

        fractpartx = modf (frequency * pt.x , &intpartx);
        fractparty = modf (frequency * pt.y , &intparty);

        MyPoint dirPt = MyPoint(intpartx, intparty) - pt;

        MyPoint triPt(cos(dirPt.x), sin(dirPt.y));
        triPt *= 0.04;

        _points[a] = tempPoints[a] + triPt;
    }*/



    ResampleCurve();
    CreateCurveVAO();
    if(_currentIter % 500 == 0)
    {
        SaveToSvg();
    }
    _currentIter++;
}


void GLWidget::PaintCurve()
{
    if(_points.size() == 0) { return; }

    int use_color_location = _shaderProgram->uniformLocation("use_color");
    _shaderProgram->setUniformValue(use_color_location, (GLfloat)1.0);

    if(SystemParams::show_points)
    {
        glPointSize(5.0f);
        _pointsVao.bind();
        glDrawArrays(GL_POINTS, 0, _points.size());
        _pointsVao.release();
    }


    //if(_drawSelPoint)
    //{
        /*
        glPointSize(15.0f);
        _selPointsVao.bind();
        glDrawArrays(GL_POINTS, 0, 1);
        _selPointsVao.release();
        */

        /*
        _rPointsVao.bind();
        glDrawArrays(GL_POINTS, 0, _rPoints.size());
        _rPointsVao.release();

        _lPointsVao.bind();
        glDrawArrays(GL_POINTS, 0, _lPoints.size());
        _lPointsVao.release();

        glLineWidth(2.0f);
        _lLinesVao.bind();
        glDrawArrays(GL_LINES, 0, _lLines.size() * 2);
        _lLinesVao.release();

        glLineWidth(2.0f);
        _rLinesVao.bind();
        glDrawArrays(GL_LINES, 0, _rLines.size() * 2);
        _rLinesVao.release();
        */
    //}

    glLineWidth(2.0f);
    _linesVao.bind();
    glDrawArrays(GL_LINES, 0, _points.size() * 2);
    _linesVao.release();
}

void GLWidget::UniformResample(std::vector<MyPoint>& oriCurve, std::vector<MyPoint>& resampleCurve, double maxDist)
{
    // fix me
    MyPoint lastPt = oriCurve[0];
    oriCurve.push_back(lastPt);

    resampleCurve.clear();

    double resample_size = maxDist;

    int iter = 0;
    double dist = 0.0;

    while(iter < oriCurve.size() - 1)
    {
        double last_dist = oriCurve[iter].Distance(oriCurve[iter + 1]);

        dist += last_dist;

        if (dist >= resample_size)
        {
            //put a point on line
            double _d = last_dist - (dist - resample_size);
            MyPoint cp(oriCurve[iter].x, oriCurve[iter].y);

            MyPoint  cp1 = MyPoint(oriCurve[iter+1].x, oriCurve[iter+1].y);

            MyPoint dirv = cp1 - cp;
            dirv = dirv * (1.0 / dirv.Length());

            MyPoint insertPt1 = cp + dirv * _d;
            resampleCurve.push_back(insertPt1);

            dist = last_dist - _d; //remaining dist

            //if remaining dist to next point needs more sampling... (within some epsilon)
            while (dist - resample_size > 1e-8 )
            {
                MyPoint insertPt2 = resampleCurve[resampleCurve.size() - 1] + dirv * resample_size;
                resampleCurve.push_back(insertPt2);
                dist -= resample_size;
            }
        }
        iter++;
    }

    //oriCurve.pop_back();
}

float GLWidget::GetRandomNumber()
{
    /*
    std::random_device rd;
    std::mt19937 e2(rd());
    std::normal_distribution<> dist(0, SystemParams::dist_std_dev);
    */
    //return ((float)std::round(_dist(_e2))) / ((float)SystemParams::dist_std_dev);
    return ((float)_dist(_e2)) / ((float)SystemParams::dist_std_dev);
}

MyPoint GLWidget::GetAttractionRepulsion(int ptIdx)
{
    std::vector<MyPoint> cPoints = _allCPoints[ptIdx];
    //GetClosestPoints(_points[ptIdx], cPoints);
    float aX = 0.0f;
    float aY = 0.0f;
    float clampVal = SystemParams::ar_clamp;
    for(size_t a = 0; a < cPoints.size(); a++)
    {
        MyPoint minVec = _points[ptIdx] - cPoints[a];
        float lVec = minVec.Length() / (SystemParams::D * GetDelta(_points[ptIdx].x, _points[ptIdx].y));
        float lj = GetLennardJones(lVec);
        MyPoint fij = (minVec / lVec ) * lj;
        float fij_length = fij.Length();

        if(fij_length > -clampVal && fij_length < clampVal)
        {
            aX += fij.x;
            aY += fij.y;
        }
    }



    return MyPoint(aX, aY) * SystemParams::f_a;
}

float GLWidget::GetLennardJones(float r)
{
    float dljr = SystemParams::sigma_l_j / r;
    //float pow_6_dljr = pow(dljr, 6.0);    // more precise
    float pow_6_dljr = fastPow(dljr, 6.0);  // approximation
    return (pow_6_dljr * pow_6_dljr) - pow_6_dljr;
}

MyPoint GLWidget::GetClosestPointToALine(MyPoint v, MyPoint w, MyPoint p)
{
    // Return minimum distance between line segment vw and point p
    double l2 = v.DistanceSquared(w);					   // i.e. |w-v|^2 -  avoid a sqrt
    if (l2 > -1e-8 && l2 < 1e-8) return v;   // v == w case

    // Consider the line extending the segment, parameterized as v + t (w - v).
    // We find projection of point p onto the line.
    // It falls where t = [(p-v) . (w-v)] / |w-v|^2
    double t = (p - v).Dot(w - v) / l2;

    if (t < 0.0)	  { return  v; }       // Beyond the 'v' end of the segment
    else if (t > 1.0) { return  w; }        // Beyond the 'w' end of the segment
    MyPoint projection = v + (w - v) * t;     // Projection falls on the segment
    return projection;
}

void GLWidget::PreparePointsVAO(std::vector<MyPoint> points, QOpenGLBuffer* ptsVbo, QOpenGLVertexArrayObject* ptsVao, QVector3D vecCol)
{
    if(ptsVao->isCreated())
    {
        ptsVao->destroy();
    }

    ptsVao->create();
    ptsVao->bind();

    QVector<VertexData> data;
    for(size_t a = 0; a < points.size(); a++)
    {
        data.append(VertexData(QVector3D(points[a].x, points[a].y,  0), QVector2D(), vecCol));
    }

    ptsVbo->create();
    ptsVbo->bind();
    ptsVbo->allocate(data.data(), data.size() * sizeof(VertexData));

    quintptr offset = 0;

    int vertexLocation = _shaderProgram->attributeLocation("vert");
    _shaderProgram->enableAttributeArray(vertexLocation);
    _shaderProgram->setAttributeBuffer(vertexLocation, GL_FLOAT, 0, 3, sizeof(VertexData));

    offset += sizeof(QVector3D);
    offset += sizeof(QVector2D);

    _shaderProgram->enableAttributeArray(_colorLocation);
    _shaderProgram->setAttributeBuffer(_colorLocation, GL_FLOAT, offset, 3, sizeof(VertexData));

    ptsVao->release();
}

void GLWidget::PrepareLinesVAO(std::vector<MyLine> lines, QOpenGLBuffer* linesVbo, QOpenGLVertexArrayObject* linesVao, QVector3D vecCol)
{
    if(linesVao->isCreated())
    {
        linesVao->destroy();
    }

    linesVao->create();
    linesVao->bind();

    QVector<VertexData> data;
    for(size_t a = 0; a < lines.size(); a++)
    {
        data.append(VertexData(QVector3D(lines[a].XA, lines[a].YA,  0), QVector2D(), vecCol));
        data.append(VertexData(QVector3D(lines[a].XB, lines[a].YB,  0), QVector2D(), vecCol));
    }

    linesVbo->create();
    linesVbo->bind();
    linesVbo->allocate(data.data(), data.size() * sizeof(VertexData));

    quintptr offset = 0;

    int vertexLocation = _shaderProgram->attributeLocation("vert");
    _shaderProgram->enableAttributeArray(vertexLocation);
    _shaderProgram->setAttributeBuffer(vertexLocation, GL_FLOAT, 0, 3, sizeof(VertexData));

    offset += sizeof(QVector3D);
    offset += sizeof(QVector2D);

    _shaderProgram->enableAttributeArray(_colorLocation);
    _shaderProgram->setAttributeBuffer(_colorLocation, GL_FLOAT, offset, 3, sizeof(VertexData));

    linesVao->release();
}

void GLWidget::GetClosestPoints()
{
    _allCPoints.clear();
    for(size_t a = 0; a < _points.size(); a++)
    {
        _allCPoints.push_back(std::vector<MyPoint>());
    }

    // KD TREE
    using namespace nanoflann;
    const float search_radius = static_cast<float>(SystemParams::kdtree_radius);
    LineCloud<float> lineCloud;
    lineCloud.lines.resize(_points.size());

    // hack
    _points.push_back(_points[0]);

    for(size_t a = 0; a < _points.size() - 1; a++)
    {
        MyPoint pt1;
        //MyPoint pt2;

        pt1 = _points[a];
        lineCloud.lines[a].x = pt1.x;
        lineCloud.lines[a].y = pt1.y;

        lineCloud.lines[a].index0 = a;
        if(a < _points.size() - 1) { /*pt2 = _points[a + 1];*/ lineCloud.lines[a].index1 = a + 1;}
        else { /*pt2 = _points[0];*/ lineCloud.lines[a].index1 = 0;}
    }

    typedef KDTreeSingleIndexAdaptor< L2_Simple_Adaptor<float, LineCloud<float> >, LineCloud<float>, 2> my_kd_tree_t;
    my_kd_tree_t   index(2 /*dim*/, lineCloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
    index.buildIndex();

    // hack
    _points.pop_back();

    for(size_t a = 0; a < _points.size(); a++)
    {
        MyPoint curPt = _points[a];
        std::vector<std::pair<size_t, float> >  ret_matches;
        nanoflann::SearchParams params;
        const float query_pt[2] = { curPt.x, curPt.y};
        const size_t nMatches = index.radiusSearch(&query_pt[0],search_radius, ret_matches, params);

        for (size_t i = 0; i < nMatches; i++)
        {
            int dataIdx = ret_matches[i].first;
            int idx1 = lineCloud.lines[dataIdx].index0;
            int idx2 = lineCloud.lines[dataIdx].index1;

            MyPoint pt1 = _points[idx1];
            MyPoint pt2 = _points[idx2];

            //MyPoint closestPt = GetClosestPointToALine(MyPoint(pt1.x, pt1.y), MyPoint(pt2.x, pt2.y), curPt);  // more precise
            MyPoint closestPt = pt1 + (pt2 - pt1) * 0.5;    // approximation
            closestPt.index = idx2;
            _allCPoints[a].push_back(closestPt);
        }
    }
}


void GLWidget::GetClosestPoints(MyPoint curPt, std::vector<MyPoint>& cPoints)
{
    // KD TREE
    using namespace nanoflann;
    const float search_radius = static_cast<float>(SystemParams::kdtree_radius);
    LineCloud<float> lineCloud;
    lineCloud.lines.resize(_points.size());

    // hack
    _points.push_back(_points[0]);

    for(size_t a = 0; a < _points.size() - 1; a++)
    {
        MyPoint pt1;
        //MyPoint pt2;

        pt1 = _points[a];
        lineCloud.lines[a].x = pt1.x;
        lineCloud.lines[a].y = pt1.y;

        lineCloud.lines[a].index0 = a;
        if(a < _points.size() - 1) { /*pt2 = _points[a + 1];*/ lineCloud.lines[a].index1 = a + 1;}
        else { /*pt2 = _points[0];*/ lineCloud.lines[a].index1 = 0;}
    }

    typedef KDTreeSingleIndexAdaptor< L2_Simple_Adaptor<float, LineCloud<float> >, LineCloud<float>, 2> my_kd_tree_t;
    my_kd_tree_t   index(2 /*dim*/, lineCloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
    index.buildIndex();

    // hack
    _points.pop_back();

    std::vector<std::pair<size_t, float> >  ret_matches;
    nanoflann::SearchParams params;
    const float query_pt[2] = { curPt.x, curPt.y};
    const size_t nMatches = index.radiusSearch(&query_pt[0],search_radius, ret_matches, params);

    for (size_t i = 0; i < nMatches; i++)
    {
        int dataIdx = ret_matches[i].first;
        int idx1 = lineCloud.lines[dataIdx].index0;
        int idx2 = lineCloud.lines[dataIdx].index1;

        MyPoint pt1 = _points[idx1];
        MyPoint pt2 = _points[idx2];

        //MyPoint closestPt = GetClosestPointToALine(MyPoint(pt1.x, pt1.y), MyPoint(pt2.x, pt2.y), curPt);  // more precise
        MyPoint closestPt = pt1 + (pt2 - pt1) * 0.5;    // approximation
        closestPt.index = idx2;
        cPoints.push_back(closestPt);
    }


}

void GLWidget::SaveToSvg()
{
    std::cout << "void GLWidget::SaveToSvg()\n";

    int wI = 1500;
    int hI = 1500;

    QSvgGenerator generator;
    generator.setFileName("image" +  QString::number(_currentIter) + ".svg");
    generator.setSize(QSize(wI, hI));
    generator.setViewBox(QRect(0, 0, wI, hI));
    generator.setTitle(tr("Bizzare Labyrinth"));
    //generator.setDescription(tr("An SVG drawing created by the SVG Generator "
    //                             "Example provided with Qt."));

    QString descriptionStr = "";
    descriptionStr += "Time (seconds): " + QString::number(_runningTime) + "\n";
    descriptionStr += "Iteration: " + QString::number(_currentIter) + "\n";
    descriptionStr += "# Points: " + QString::number(_points.size()) + "\n";
    descriptionStr += "f_b: " + QString::number(SystemParams::f_b) + "\n";
    descriptionStr += "f_f: " + QString::number(SystemParams::f_f) + "\n";
    descriptionStr += "f_a: " + QString::number(SystemParams::f_a) + "\n";
    descriptionStr += "ar_clamp: " + QString::number(SystemParams::ar_clamp) + "\n";
    descriptionStr += "D: " + QString::number(SystemParams::D) + "\n";
    descriptionStr += "sigma_l_j: " + QString::number(SystemParams::sigma_l_j) + "\n";
    descriptionStr += "search_radius: " + QString::number(SystemParams::search_radius) + "\n";

    generator.setDescription(descriptionStr);

    QPainter painter;
    painter.begin(&generator);

    // draw
    painter.setClipRect(QRect(0, 0, _img_width, _img_height));
    painter.setPen(QPen(Qt::black, 1.0));

    int mulFactor = 10;
    int xFactor = (wI / 2) - (this->_img_width / 2 * mulFactor);
    int yFactor = (hI / 2) - (this->_img_width / 2 * mulFactor);

    for(size_t a = 0; a < _points.size(); a++)
    {
        MyPoint startPt = _points[a];
        MyPoint endPt;

        if(a == _points.size() - 1)
        {
            endPt = _points[0];
        }
        else
        {
            endPt = _points[a + 1];
        }

        startPt *= mulFactor;
        endPt *= mulFactor;

        startPt.x += xFactor;
        startPt.y += yFactor;
        endPt.x += xFactor;
        endPt.y += yFactor;

        painter.drawLine(startPt.x, startPt.y, endPt.x, endPt.y);
    }

    painter.end();
}

void GLWidget::PrepareImageVAO(QOpenGLTexture* tex, QOpenGLBuffer* vbo, QOpenGLVertexArrayObject* vao)
{
    // ~~~ create vao for the input image ~~~
    vao->create();
    vao->bind();

    QVector<VertexData> imageVertices;
    imageVertices.append(VertexData(QVector3D(0.0,        0.0,          0.0f), QVector2D(0, 0)));
    imageVertices.append(VertexData(QVector3D(_img_width, 0.0,          0.0f), QVector2D(1, 0)));
    imageVertices.append(VertexData(QVector3D(_img_width, _img_height,  0.0f), QVector2D(1, 1)));
    imageVertices.append(VertexData(QVector3D(0.0,        _img_height,  0.0f), QVector2D(0, 1)));

    vbo->create();
    vbo->bind();
    vbo->allocate(imageVertices.data(), 4 * sizeof(VertexData));

    //tex->setMinificationFilter(QOpenGLTexture::Nearest);
    //tex->setMagnificationFilter(QOpenGLTexture::Linear);
    //_shaderProgram->setAttributeValue("base_texture", tex->textureId());

    // Offset for position
    quintptr offset = 0;

    // vertex
    int vertexLocation = _shaderProgram->attributeLocation("vert");
    _shaderProgram->enableAttributeArray(vertexLocation);
    _shaderProgram->setAttributeBuffer(vertexLocation, GL_FLOAT, offset, 3, sizeof(VertexData));

    offset += sizeof(QVector3D);

    // uv
    int texcoordLocation = _shaderProgram->attributeLocation("uv");
    _shaderProgram->enableAttributeArray(texcoordLocation);
    _shaderProgram->setAttributeBuffer(texcoordLocation, GL_FLOAT, offset, 2, sizeof(VertexData));

    vao->release();
}


void GLWidget::SetImage(QString img)
{
    _imgColor.load(img);
    _imgOriginal = LoadImageAsGrayscale(img);

    // size
    this->_img_width = _imgOriginal.width() ;
    this->_img_height =  _imgOriginal.height() ;

    std::cout << "image loaded. width: " << _img_width << ", height: " << _img_height << "\n";

    _imageTexture = new QOpenGLTexture(_imgOriginal);

    PrepareImageVAO(_imageTexture, &_imageVbo, &_imageVao);

    CalculateGradient();

}

void GLWidget::CalculateGradient()
{
    // calculate gradient
    _imgGradientX = QImage(_img_width, _img_height, QImage::Format_ARGB32);
    _imgGradientY = QImage(_img_width, _img_height, QImage::Format_ARGB32);
    _imgMagnitude = QImage(_img_width, _img_height, QImage::Format_ARGB32);
    for (int a = 0; a < _img_height; a++)
    {
        uchar* scanX = _imgGradientX.scanLine(a);
        uchar* scanY = _imgGradientY.scanLine(a);
        uchar* scanM = _imgMagnitude.scanLine(a);
        int depth = 4;
        for (int b = 0; b < _img_width; b++)
        {
            QRgb* rgbpixelX = reinterpret_cast<QRgb*>(scanX + b * depth);
            *rgbpixelX = QColor(127, 127, 127).rgba();   // gray 0.5

            QRgb* rgbpixelY = reinterpret_cast<QRgb*>(scanY + b * depth);
            *rgbpixelY = QColor(127, 127, 127).rgba();   // gray 0.5

            QRgb* rgbpixelM = reinterpret_cast<QRgb*>(scanM + b * depth);
            *rgbpixelM = QColor(127, 127, 127).rgba();   // gray 0.5

        }
    }

    for(size_t a = 0; a < _img_width; a++)
    {
        _gradientX.push_back(std::vector<float>(_img_height));
        _gradientY.push_back(std::vector<float>(_img_height));
        _magnitude.push_back(std::vector<float>(_img_height));

        for(size_t b = 0; b < _img_height; b++)
           {
            _gradientX[a][b] = 0.5f;
            _gradientY[a][b] = 0.5f;
            _magnitude[a][b] = 0.5f;
        }
    }

    float kernelx[3][3] = {{-1, -2, -1},
                           {0,  0,  0},
                           {1,  2,  1}};


    float kernely[3][3] = {{-1, 0, 1},
                           {-2, 0, 2},
                           {-1, 0, 1}};

    for(size_t a = 1; a < _img_width - 1; a++)
    {
        for(size_t b = 1; b < _img_height - 1; b++)
        {
            float gx = 0.0;
            float gy = 0.0;

            for(size_t i = 0; i < 3; i++)
            {
                for(size_t j = 0; j < 3; j++)
                {
                    int an = a + i - 1;
                    int bn = b + j - 1;

                    float pixVal = ((float)qGray(_imgOriginal.pixel(an, bn))) / 255.0;
                    gx += pixVal * kernelx[i][j];
                    gy += pixVal * kernely[i][j];
                }
             }


            float mag = sqrt(gx * gx + gy * gy);
            _gradientX[a][b] = gx;
            _gradientY[a][b] = gy;
            _magnitude[a][b] = mag;

            int gx255 = (int)(gx * 255.0);
            _imgGradientX.setPixel(a, b, qRgb(gx255, gx255, gx255));
            int gy255 = (int)(gy * 255.0);
            _imgGradientY.setPixel(a, b, qRgb(gy255, gy255, gy255));
            int mag255 = (int)(mag * 255.0);
            _imgMagnitude.setPixel(a, b, qRgb(mag255, mag255, mag255));
        }
    }

    _magnitudeTexture = new QOpenGLTexture(_imgMagnitude);
    PrepareImageVAO(_magnitudeTexture, &_magnitudeVbo, &_magnitudeVao);
}

QImage GLWidget::LoadImageAsGrayscale(QString img)
{
    QImage oriImage;
    oriImage.load(img);
    for (int i = 0; i < oriImage.height(); i++)
    {
        uchar* scan = oriImage.scanLine(i);
        int depth = 4;
        for (int jj = 0; jj < oriImage.width(); jj++) {

            QRgb* rgbpixel = reinterpret_cast<QRgb*>(scan + jj * depth);
            int gray = qGray(*rgbpixel);
            *rgbpixel = QColor(gray, gray, gray).rgba();
        }
    }
    return oriImage;
}
