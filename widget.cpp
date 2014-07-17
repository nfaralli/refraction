#include "widget.h"
#include <QPainter>
#include <QMouseEvent>
#include <QWheelEvent>
#include <QFileDialog>
#include <QApplication>
#include <QDesktopWidget>

Widget::Widget(QWidget *parent) : QWidget(parent) {

  float alpha=0.7;
  int   xScreen=QApplication::desktop()->width();
  int   yScreen=QApplication::desktop()->height();
  //biggest aspect ratio=(4/3)^3. if bigger, most likely a dual screen
  if(((float)xScreen)/yScreen > (64./27.*1.1))
    xScreen/=2;
  if(((float)xScreen)/yScreen < (4./3./1.1))
    yScreen/=2;
  resize((int)(xScreen*alpha),(int)(yScreen*alpha));
  move((int)(xScreen*(1-alpha)/2),(int)(yScreen*(1-alpha)/2));
  setWindowTitle("refraction");
  setBackgroundRole(QPalette::WindowText);
  setAutoFillBackground(true);
  m_bg_color = this->palette().color(this->backgroundRole()).rgb();
  m_pressed_button = Qt::NoButton;
  m_shift_key_pressed = false;
  m_T_key_pressed = false;

  int radius=100;
  m_filter_method=getSimpleLensFilter;
  m_lens_params.n1=1.0;
  m_lens_params.n2=1.52; //refractive index of crown glass
  m_lens_params.radius=radius;
  m_lens_params.thickness=20;
  m_lens_params.height=0;
  m_lens_params.lens_type=ONE_CURVED_SIDE_TOP;
  m_filter=m_filter_method(&m_lens_params);

  m_image = QImage("fluffinette.jpg");
  if(m_image.isNull()){
    loadImage();
  }
  m_image_offset = QPoint(0, 0);
  m_image_sphere = QImage(2*radius+1, 2*radius+1, QImage::Format_ARGB32);
  m_image_sphere_offset = QPoint(this->size().width()/2 - radius, this->size().height()/2 - radius);
  m_offset_delta = m_image_sphere_offset - m_image_offset;
}

void Widget::loadImage(){
  QString fileName = QFileDialog::getOpenFileName(
      this,QString(),QString(),QString("jpeg (*.jpg)"));
  m_image = QImage(fileName);
  if(m_image.isNull()){
    loadImage();
  }
  else{
    update();
  }
}

QRgb Widget::getColor(int i, int j) {
  int index, k, num_pts;
  int ii, jj;
  float r, g, b;
  Point *pts, *pt;
  float *coefs, *coef;
  QRgb color;


  // should check range of i and j first
  index=i+m_filter->dim_x*j;
  num_pts=m_filter->num_pts[index];
  if(num_pts==0){
    ii=i+m_offset_delta.x();
    jj=j+m_offset_delta.y();
    if(ii>=0 && ii<m_image.width() && jj>=0 && jj<m_image.height())
      return m_image.pixel(ii, jj);
    else
      return m_bg_color;
  }
  pts=m_filter->pts[index];
  coefs=m_filter->coefs[index];
  r=g=b=0;
  for(pt=pts, coef=coefs, k=0; k<num_pts; k++, pt++, coef++){
    ii=pt->x+m_offset_delta.x()+m_lens_params.radius;
    jj=pt->y+m_offset_delta.y()+m_lens_params.radius;
    if(ii>=0 && ii<m_image.width() && jj>=0 && jj<m_image.height())
      color=m_image.pixel(ii, jj);
    else
      color=m_bg_color;
    r+=qRed(color)*(*coef);
    g+=qGreen(color)*(*coef);
    b+=qBlue(color)*(*coef);
  }
  return qRgb((int)r, (int)g, (int)b);
}

void Widget::setSphere(){
  int i, j;
  for(j=0; j<m_filter->dim_y; j++){
    for(i=0; i<m_filter->dim_x; i++){
      m_image_sphere.setPixel(i, j, getColor(i, j));
    }
  }
}

void Widget::paintEvent(QPaintEvent * /* event */){
  QPainter painter(this);
  setSphere();
  painter.drawImage(m_image_offset, m_image);
  painter.drawImage(m_image_sphere_offset, m_image_sphere);
}

void Widget::mousePressEvent(QMouseEvent *e){
  switch(e->button()){
  case Qt::LeftButton:
    m_left_pos = e->pos();
    m_image_sphere_offset = m_left_pos-QPoint(m_lens_params.radius, m_lens_params.radius);
    m_offset_delta = m_image_sphere_offset - m_image_offset;
    break;
  case Qt::RightButton:
    m_right_pos = e->pos();
    break;
  default:
      break;
  }
  m_pressed_button = e->button();
  update();
}

void Widget::mouseReleaseEvent(QMouseEvent * /* event */){
  m_pressed_button = Qt::NoButton;
}

void Widget::mouseMoveEvent(QMouseEvent *e){
  switch(m_pressed_button){
  case Qt::RightButton:
    m_image_offset += e->pos() - m_right_pos;
    m_image_sphere_offset += e->pos() - m_right_pos;
    m_right_pos = e->pos();
    break;
  case Qt::LeftButton:
    m_image_sphere_offset += e->pos() - m_left_pos;
    m_offset_delta += e->pos() - m_left_pos;
    m_left_pos = e->pos();
    break;
  default:
    break;
  }
  update();
}

void Widget::wheelEvent(QWheelEvent *e){
  if (e->orientation() == Qt::Horizontal){
    m_image_offset.rx() -= e->delta();
    m_image_sphere_offset.rx() -= e->delta();
  }
  else{
    m_image_offset.ry() += e->delta();
    m_image_sphere_offset.ry() += e->delta();
  }
  update();
}

void Widget::keyPressEvent(QKeyEvent *e){
  bool update_filter=false;
  QPoint offset;

  switch(e->key()){
  case Qt::Key_Plus:
    if(m_shift_key_pressed){
      m_lens_params.height+=5;
    }
    else if(m_T_key_pressed){
      if(m_lens_params.thickness+1<=m_lens_params.radius){
        m_lens_params.thickness++;
      }
    }
    else{
      m_lens_params.radius+=5;
      offset=-QPoint(5, 5);
    }
    update_filter=true;
    break;
  case Qt::Key_Minus:
    if(m_shift_key_pressed){
      m_lens_params.height-=5;
    }
    else if(m_T_key_pressed){
      if(m_lens_params.thickness>=2){
        m_lens_params.thickness--;
      }
    }
    else{
      if(m_lens_params.radius>=10 && m_lens_params.thickness<=(m_lens_params.radius-5)){
        m_lens_params.radius-=5;
        offset=QPoint(5, 5);
      }
    }
    update_filter=true;
    break;
  case Qt::Key_Shift:
    m_shift_key_pressed=true;
    break;
  case Qt::Key_L:
    m_lens_params.lens_type=(LensType)((m_lens_params.lens_type+1)%NUM_LENS_TYPE);
    update_filter=true;
    break;
  case Qt::Key_S:
    if(m_filter_method == getLensFilter)
      m_filter_method=getSimpleLensFilter;
    else
      m_filter_method=getLensFilter;
    update_filter=true;
    break;
  case Qt::Key_C:
    setCursor(Qt::BlankCursor);
    break;
  case Qt::Key_O:
    loadImage();
    break;
  case Qt::Key_T:
    m_T_key_pressed=true;
    break;
  }
  if(update_filter){
    int radius = m_lens_params.radius;
    freeFilter(m_filter);
    m_filter=m_filter_method(&m_lens_params);
    m_image_sphere = QImage(2*radius+1, 2*radius+1, QImage::Format_ARGB32);
    m_image_sphere_offset+=offset;
    m_offset_delta+=offset;
    update();
  }
}

void Widget::keyReleaseEvent(QKeyEvent *e){
  switch(e->key()){
  case Qt::Key_Shift:
     m_shift_key_pressed=false;
    break;
  case Qt::Key_C:
    setCursor(Qt::ArrowCursor);
    break;
  case Qt::Key_T:
     m_T_key_pressed=false;
    break;
  }
}
