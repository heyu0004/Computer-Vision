// generated by Fast Light User Interface Designer (fluid) version 1.0011

#ifndef HelpPageUI_h
#define HelpPageUI_h
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Output.H>
#include <FL/Fl_Button.H>

class HelpPageUI {
public:
  HelpPageUI();
  Fl_Window *mainWindow;
  Fl_Output *helpText;
private:
  inline void cb_OK_i(Fl_Button*, void*);
  static void cb_OK(Fl_Button*, void*);
public:
  void show();
  void hide();
};
#endif
