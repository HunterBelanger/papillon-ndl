#include <PapillonNDL/polynomial_1d.hpp>
#include <gtest/gtest.h>

namespace {
  using namespace pndl;

  TEST(Polynomial1D, Order) {
    std::vector<double> coeffs {3.,4.,5.,6.};
    Polynomial1D poly(coeffs);

    EXPECT_EQ(poly.order(), coeffs.size() - 1);

    std::vector<double> coeffs2 {3.,4.,5.,6.,2.,1.,3.5,6.5};
    Polynomial1D poly2(coeffs2);

    EXPECT_EQ(poly2.order(), coeffs2.size() - 1);
  }

  TEST(Polynomial1D, Coefficients) {
    std::vector<double> coeffs {1.1,2.2,3.3,4.4};
    Polynomial1D poly(coeffs);

    for(size_t i = 0; i <= poly.order(); i++) {
      EXPECT_DOUBLE_EQ(coeffs[i], poly.coefficient(i)); 
    }
  }

  TEST(Polynomial1D, Evaluation) {
    std::vector<double> coeffs {1.1,2.2,3.3,4.4};
    Polynomial1D poly(coeffs);

    double x = 0.;
    double y = coeffs[0];
    EXPECT_DOUBLE_EQ(y, poly(x));

    x = 1.;
    y = coeffs[0] + coeffs[1] + coeffs[2] + coeffs[3];
    EXPECT_DOUBLE_EQ(y, poly(x));

    x = 2.;
    y = 53.900000000000006;
    EXPECT_DOUBLE_EQ(y, poly(x));

    x = 5.;
    y = 644.6;
    EXPECT_DOUBLE_EQ(y, poly(x));

    x = 20.;
    y = 36565.1;
    EXPECT_DOUBLE_EQ(y, poly(x));
  }

  TEST(Polynomial1D, Integration) {
    std::vector<double> coeffs {1.1,2.2,3.3,4.4};
    Polynomial1D poly(coeffs);

    double x_low = 1.;
    double x_hi = 5.; 
    double i = 853.6;
    EXPECT_DOUBLE_EQ(i, poly.integrate(x_low, x_hi));

    x_low = -7.8;
    x_hi = 22.7;
    i = 301926.74985;
    EXPECT_DOUBLE_EQ(i, poly.integrate(x_low, x_hi));

    x_hi = -7.8;
    x_low = 22.7;
    EXPECT_DOUBLE_EQ(-i, poly.integrate(x_low, x_hi));
  }

}
