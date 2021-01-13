#include <PapillonNDL/interpolation.hpp>
#include <gtest/gtest.h>

namespace {
  using namespace pndl;

  TEST(Interpolation, Histogram) {
    Histogram interp;

    double x1 = 8.;
    double y1 = 4.;

    double x2 = 10.;
    double y2 = 6.;

    EXPECT_DOUBLE_EQ(y1, interp.interpolate(x1, x1, y1, x2, y2));
    
    double x = 9.;
    double y = y1;
    EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));

    x = x2 - 0.000000001;
    EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));

    EXPECT_DOUBLE_EQ(y1, interp.interpolate(x2, x1, y1, x2, y2));
  }

  TEST(Interpolation, LinearLinear) {
    LinLin interp;
    // Here, y is linear in x
    // y = ((x - x1)/(x2 - x1))*(y2 - y1) + y1

    double x1 = 8.;
    double y1 = 0.;

    double x2 = 10.;
    double y2 = 6.;
    EXPECT_DOUBLE_EQ(y1, interp.interpolate(x1, x1, y1, x2, y2));
    
    double x = 9.;
    double y = 3.;
    EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));
    
    EXPECT_DOUBLE_EQ(y2, interp.interpolate(x2, x1, y1, x2, y2));
  }

  TEST(Interpolation, LinearLog) {
    // Here, y is linear in log(x)
    // y = (log(x/x1)/log(x2/x1))*(y2 - y1) + y1
    LinLog interp;
    
    double x1 = 1.;
    double y1 = 1.;

    double x2 = 4.;
    double y2 = 3.;
    EXPECT_DOUBLE_EQ(y1, interp.interpolate(x1, x1, y1, x2, y2));
    
    double x = 2.;
    double y = 2.;
    EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));

    x = 3.;
    y = 2.584962500721156;
    EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));
    
    EXPECT_DOUBLE_EQ(y2, interp.interpolate(x2, x1, y1, x2, y2));
  }

  TEST(Interpolation, LogLinear) {
    // Here, y is linear in log(x)
    // log(y) = ((x - x1)/(x2 - x1))*log(y2/y1) + log(y1)
    LogLin interp;
    
    double x1 = 1.;
    double y1 = 1.;

    double x2 = 4.;
    double y2 = 3.;
    EXPECT_DOUBLE_EQ(y1, interp.interpolate(x1, x1, y1, x2, y2));

    double x = 2.;
    double y = 1.4422495703074083;
    EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));

    x = 3.;
    y = 2.080083823051904;
    EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));
    
    
    EXPECT_DOUBLE_EQ(y2, interp.interpolate(x2, x1, y1, x2, y2));
  }

  TEST(Interpolation, LogLog) {
    // Here, log(y) is linear in log(x)
    // log(y) = (log(x/x1)/log(x2/x1))*log(y2/y1) + log(y1)
    LogLog interp;
    
    double x1 = 1.;
    double y1 = 1.;

    double x2 = 4.;
    double y2 = 3.;
    EXPECT_DOUBLE_EQ(y1, interp.interpolate(x1, x1, y1, x2, y2));

    double x = 2.;
    double y =1.7320508075688774 ;
    EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));

    x = 3.;
    y = 2.388414221757005;
    EXPECT_DOUBLE_EQ(y, interp.interpolate(x, x1, y1, x2, y2));
    
    
    EXPECT_DOUBLE_EQ(y2, interp.interpolate(x2, x1, y1, x2, y2));
  }

}
