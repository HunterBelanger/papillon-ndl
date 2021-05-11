#include <PapillonNDL/pctable.hpp>
#include <gtest/gtest.h>
#include <vector>

namespace pndl {
namespace {

TEST(PCTable, Construction) {
  /// v1 unsorted
  std::vector<double> v1 {1., 3., 2.};
  std::vector<double> p1 {0.7, 0.3, 0.3};
  std::vector<double> c1 {0., 0.7, 1.};
  EXPECT_THROW(PCTable(v1, p1, c1, Interpolation::LinLin), PNDLException);

  // pdf has negative values
  std::vector<double> v2 {1., 2., 3.};
  std::vector<double> p2 {0.7, -0.3, 0.3};
  std::vector<double> c2 {0., 0.7, 1.};
  EXPECT_THROW(PCTable(v2, p2, c2, Interpolation::LinLin), PNDLException);

  // cdf not monotonic 
  std::vector<double> v3 {1., 2., 3., 4.};
  std::vector<double> p3 {0.7, 0.3, 0.3, 0.6};
  std::vector<double> c3 {0., 0.3, 0.2, 1.};
  EXPECT_THROW(PCTable(v3, p3, c3, Interpolation::LinLin), PNDLException);

  // cdf doesnt end at 1
  std::vector<double> v4 {1., 2., 3., 4.};
  std::vector<double> p4 {0.7, 0.3, 0.3, 0.6};
  std::vector<double> c4 {0., 0.3, 0.2, 0.9};
  EXPECT_THROW(PCTable(v4, p4, c4, Interpolation::LinLin), PNDLException);

  // Good grids, bad interpolation
  std::vector<double> v5 {1., 2., 3., 4.};
  std::vector<double> p5 {0.7, 0.3, 0.3, 0.6};
  std::vector<double> c5 {0., 0.3, 0.6, 1.};
  EXPECT_THROW(PCTable(v5, p5, c5, Interpolation::LinLog), PNDLException);
  EXPECT_THROW(PCTable(v5, p5, c5, Interpolation::LogLin), PNDLException);
  EXPECT_THROW(PCTable(v5, p5, c5, Interpolation::LogLog), PNDLException);

  // Good everything
  EXPECT_NO_THROW(PCTable(v5, p5, c5, Interpolation::Histogram));
  EXPECT_NO_THROW(PCTable(v5, p5, c5, Interpolation::LinLin));
}

TEST(PCTable, SampleValue) {
  // Histogram
  std::vector<double> vh {1., 2., 3.};
  std::vector<double> ph {0.7, 0.3, 0.3};
  std::vector<double> ch {0., 0.7, 1.};
  PCTable hist(vh, ph, ch, Interpolation::Histogram);

  EXPECT_DOUBLE_EQ(hist.sample_value(0.7), 2.);
  EXPECT_DOUBLE_EQ(hist.sample_value(0.5), 1. + 5./7.);
  EXPECT_DOUBLE_EQ(hist.sample_value(0.8), 2. + 1./3.);
  EXPECT_DOUBLE_EQ(hist.sample_value(1.), 3.);

  // Linear
  std::vector<double> vl {1., 2., 3., 4.};
  std::vector<double> pl {0., 0.25, 0.25, 1.};
  std::vector<double> cl {0., 0.125, 0.375, 1.};
  PCTable lin(vl, pl, cl, Interpolation::LinLin);

  EXPECT_DOUBLE_EQ(lin.sample_value(0.125), 2.);
  EXPECT_DOUBLE_EQ(lin.sample_value(0.375), 3.);
  EXPECT_DOUBLE_EQ(lin.sample_value(0.03125), 1.5);
  EXPECT_DOUBLE_EQ(lin.sample_value(0.2), 2.3);
}

TEST(PCTable, PDFEvaluation) {
  // Histogram
  std::vector<double> vh {1., 2., 3.};
  std::vector<double> ph {0.7, 0.3, 0.3};
  std::vector<double> ch {0., 0.7, 1.};
  PCTable hist(vh, ph, ch, Interpolation::Histogram);

  EXPECT_DOUBLE_EQ(hist.pdf(1.), 0.7);
  EXPECT_DOUBLE_EQ(hist.pdf(2.), 0.3);
  EXPECT_DOUBLE_EQ(hist.pdf(1.5), 0.7);
  EXPECT_DOUBLE_EQ(hist.pdf(2.9), 0.3);

  // Linear
  std::vector<double> vl {1., 2., 3., 4.};
  std::vector<double> pl {0., 0.25, 0.25, 1.};
  std::vector<double> cl {0., 0.125, 0.375, 1.};
  PCTable lin(vl, pl, cl, Interpolation::LinLin);

  EXPECT_DOUBLE_EQ(lin.pdf(1.5), 0.125);
  EXPECT_DOUBLE_EQ(lin.pdf(2.), 0.25);
  EXPECT_DOUBLE_EQ(lin.pdf(2.5), 0.25);
  EXPECT_DOUBLE_EQ(lin.pdf(3.75), 0.8125);
}

TEST(PCTable, MinMaxValue) {
  std::vector<double> vl {-2.45, 2., 3., 48.};
  std::vector<double> pl {0., 0.25, 0.25, 1.};
  std::vector<double> cl {0., 0.125, 0.375, 1.};
  PCTable lin(vl, pl, cl, Interpolation::LinLin);

  EXPECT_DOUBLE_EQ(lin.max_value(), 48.);
  EXPECT_DOUBLE_EQ(lin.min_value(), -2.45);
}

TEST(PCTable, Size) {
  std::vector<double> vl {-2.45, 2., 3., 48.};
  std::vector<double> pl {0., 0.25, 0.25, 1.};
  std::vector<double> cl {0., 0.125, 0.375, 1.};
  PCTable lin(vl, pl, cl, Interpolation::LinLin);

  EXPECT_EQ(lin.size(), vl.size());
}

TEST(PCTable, ValuesGrid) {
  std::vector<double> vl {-2.45, 2., 3., 48.};
  std::vector<double> pl {0., 0.25, 0.25, 1.};
  std::vector<double> cl {0., 0.125, 0.375, 1.};
  PCTable lin(vl, pl, cl, Interpolation::LinLin);

  const auto& vals = lin.values();

  for (std::size_t i = 0; i < vals.size(); i++) {
    EXPECT_DOUBLE_EQ(vals[i], vl[i]);
  }
}

TEST(PCTable, PDFGrid) {
  std::vector<double> vl {-2.45, 2., 3., 48.};
  std::vector<double> pl {0., 0.25, 0.25, 1.};
  std::vector<double> cl {0., 0.125, 0.375, 1.};
  PCTable lin(vl, pl, cl, Interpolation::LinLin);

  const auto& vals = lin.pdf();

  for (std::size_t i = 0; i < vals.size(); i++) {
    EXPECT_DOUBLE_EQ(vals[i], pl[i]);
  }
}

TEST(PCTable, CDFGrid) {
  std::vector<double> vl {-2.45, 2., 3., 48.};
  std::vector<double> pl {0., 0.25, 0.25, 1.};
  std::vector<double> cl {0., 0.125, 0.375, 1.};
  PCTable lin(vl, pl, cl, Interpolation::LinLin);

  const auto& vals = lin.cdf();

  for (std::size_t i = 0; i < vals.size(); i++) {
    EXPECT_DOUBLE_EQ(vals[i], cl[i]);
  }
}

TEST(PCTable, Interpolation) {
  std::vector<double> vl {-2.45, 2., 3., 48.};
  std::vector<double> pl {0., 0.25, 0.25, 1.};
  std::vector<double> cl {0., 0.125, 0.375, 1.};
  PCTable lin(vl, pl, cl, Interpolation::LinLin);
  EXPECT_EQ(Interpolation::LinLin, lin.interpolation());

  PCTable hist(vl, pl, cl, Interpolation::Histogram);
  EXPECT_EQ(Interpolation::Histogram, hist.interpolation());
}

}
}