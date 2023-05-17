#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <float.h>
#include <utility>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_pushButton_clicked();

private:
    Ui::MainWindow *ui;

private:
    using vec = std::vector<double>;
    using vecvec = std::vector<std::vector<double>>;

    bool idx_is_part_of_area(size_t i, size_t j);
    bool coord_is_part_of_area(double x, double y);

    template <typename func>
    std::pair<double, double> get_minmax(double a, double b, double c, double d, double xstep, double ystep, func f);

    void init_v(vecvec& v, double h, double k);

    double u_test(double i, double j, size_t _n, size_t _m) {
        const double h = (b - a) / _n;
        const double k = (d - c) / _m;
        double x = a + i * h;
        double y = c + j * k;
        return std::exp(1 - x * x - y * y);
    }

    double mu1_test(double y) {
        return  std::exp(-y * y); //x = a = -1
    }

    double mu2_test(double x) {
        return std::exp(-x * x); //y = d = 1
    }

    double mu3_test(double y) {
        return std::exp(1 - y * y); //x = (b-a)/2 = 0
    }

    double mu4_test(double x) {
        return std::exp(1 -x * x); //y = (d-c)/2 = 0
    }

    double mu5_test(double y) {
        return std::exp(-y * y); //x = b = 1
    }

    double mu6_test(double x) {
        return std::exp(-x * x); //y = c = -1
    }

    double f_test(double x, double y) {
        return (-4) * (1 - x * x - y * y) * std::exp(1 - x * x - y * y);
    }

private:
    size_t Nmax;
    double eps;
    size_t n, m;
    const double a = -1.0;
    const double b = 1.0;
    const double c = -1.0;
    const double d = 1.0;
};

template <typename func>
std::pair<double, double> MainWindow::get_minmax(double a, double b, double c, double d, double xstep, double ystep, func get_val)
{
    double max = DBL_MAX, min = max;

    for (double x = a; x <= b; x += xstep)
    {
        for (double y = c; y <= d; y += ystep)
        {
            if (coord_is_part_of_area(x, d))
            {
                size_t i = size_t((x - a) / xstep);
                size_t j = size_t((y - c) / ystep);

                double val = get_val(i, j, x, y);
                if (val > max) max = val;
                else if (val < min) min = val;
            }
        }
    }

    return {min, max};
}

#endif // MAINWINDOW_H
