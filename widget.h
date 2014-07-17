#ifndef __WIDGET_H__
#define __WIDGET_H__

#include <QWidget>
#include "utils.h"

class Widget : public QWidget {

public:
  Widget(QWidget *parent=0);

protected:
  void paintEvent(QPaintEvent *);
  void mousePressEvent(QMouseEvent *);
  void mouseMoveEvent(QMouseEvent *);
  void mouseReleaseEvent(QMouseEvent *);
  void wheelEvent(QWheelEvent *);
  void keyPressEvent(QKeyEvent *);
  void keyReleaseEvent(QKeyEvent *);

private:
  void loadImage();
  QRgb getColor(int, int);
  void setSphere();

  Filter*(*m_filter_method)(LensParams*);
  LensParams m_lens_params;
  Filter *m_filter;
  QImage m_image;
  QImage m_image_sphere;
  QPoint m_image_offset;
  QPoint m_image_sphere_offset;
  QPoint m_offset_delta;
  QPoint m_left_pos;
  QPoint m_right_pos;
  Qt::MouseButton m_pressed_button;
  QRgb m_bg_color;
  bool m_shift_key_pressed;
  bool m_T_key_pressed;
};

#endif // __WIDGET_H__
