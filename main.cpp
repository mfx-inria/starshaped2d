#include <vector>
#include <queue>
#include <limits>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

template<typename T>
class Vec2 {
public:
    Vec2(T xx) : x(xx), y(xx) {}
    Vec2(T xx, T yy) : x(xx), y(yy) {}
    Vec2 operator + (const Vec2 &v) const { return Vec2(x + v.x, y + v.y); }
    Vec2 operator - (const Vec2 &v) const { return Vec2(x - v.x, y - v.y); }
    Vec2 operator * (const T &r) const { return Vec2(x * r, y * r); }
    Vec2 operator / (const T &r) const { return Vec2(x / r, y / r); }
    T dotProduct(const Vec2<T> &v) const { return x * v.x + y * v.y; }
    Vec2& operator /= (const T &r) { x /= r, y /= r; return *this; }
    T length() const{ return sqrt(x * x + y * y); }
    T& operator [] (int i) { return (&x)[i]; }
    friend std::ostream& operator << (std::ostream &s, const Vec2<T> &v) { return s << '[' << v.x << ' ' << v.y << ']';}
    T x, y;
};
typedef Vec2<double> Vec2d;
const int neigh_indices[] = {0, 1, 0, -1, 1, 0, -1, 0}; // 4-connected neighbourhood indices
const double epsilon = std::numeric_limits<double>::epsilon();
enum LatticeType{DiagonalSquare, Triangular, Honeycomb};
typedef std::pair<double, std::pair<unsigned int, unsigned int> > GrowthInfo;

class CompareVoronoiGrowth {
public:
    bool operator()(GrowthInfo n1, GrowthInfo n2) const {
        if (n1.first > n2.first) {
            return true;
        } else {
            if (n1.first < n2.first) return false;
            return (n1.second.first < n2.second.first); // lexicographic comparison when at equal distance
        }
    }
};

void line_through_two_points(const Vec2d& p, const Vec2d& q, double& a, double& b, double& c){
    a = q.y - p.y;
    b = (p.x - q.x);
    c = -(a * p.x + b * p.y);
}

double star_shaped_polygonal_length(const std::vector<Vec2d>& star_shaped_polygon, const Vec2d& ray_direction) {
    /** Return the length from a ray centered at the origin (ray_direction * t), and the set of straight line segments
    defining the star-shaped polygon.*/
    if (ray_direction.x == 0.0 && ray_direction.y == 0.0) return 0.0;
    double min_length = std::numeric_limits<double>::max();
    for (unsigned int i = 0; i < star_shaped_polygon.size(); i++) {
        Vec2d p1 = star_shaped_polygon[i];
        Vec2d p2 = star_shaped_polygon[(i + 1) % star_shaped_polygon.size()];
        double la, lb, lc;
        line_through_two_points(p1, p2, la, lb, lc); // line equation with form ax + by + c = 0
        double t = (-lc / (Vec2d(la, lb).dotProduct(ray_direction))); // intersection point of ray and line
        if (t >= 0) {
            Vec2d inter_point = ray_direction * t;
            double length = inter_point.length();
            if (length < min_length) {
                // check that intersection point is in segment (needed due to possible non convexity of the polygon)
                if (!(inter_point[0] >= std::min(p1[0], p2[0]) - epsilon)) continue;
                if (!(inter_point[0] <= std::max(p1[0], p2[0]) + epsilon)) continue;
                if (!(inter_point[1] >= std::min(p1[1], p2[1]) - epsilon)) continue;
                if (!(inter_point[1] <= std::max(p1[1], p2[1]) + epsilon)) continue;
                min_length = length;
            }
        }
    }
    return min_length;
}

double get_distance(Vec2d& coordinate, Vec2d& site, std::vector<Vec2d>& star_shaped_polygon) {
    // return (coordinate - site).length(); // Euclidean distance
    return (coordinate - site).length() / star_shaped_polygonal_length(star_shaped_polygon, coordinate - site); // star-shaped distance
}

Vec2d get_coordinate(unsigned int x, unsigned int y, unsigned int image_size) {
    Vec2d coordinate(static_cast<double>(x) + 0.5, static_cast<double>(y) + 0.5);
    coordinate /= static_cast<double>(image_size);
    return coordinate;
}

void compute_voronoi_growth_labeling(unsigned int image_size, std::vector<Vec2d>& sites, std::vector<unsigned int>& image_closest_site, int max_growth_length, std::vector<Vec2d>& star_shaped_polygon) {
    std::priority_queue<GrowthInfo, std::vector<GrowthInfo >, CompareVoronoiGrowth> queue;
    for (unsigned int s = 0; s < sites.size(); s++){
        Vec2d grid_pos = sites[s] * static_cast<double>(image_size);
        unsigned int x = static_cast<unsigned int>(std::floor(grid_pos.x));
        unsigned int y = static_cast<unsigned int>(std::floor(grid_pos.y));
        Vec2d coordinate = get_coordinate(x, y, image_size);
        queue.push(std::make_pair(get_distance(coordinate, sites[s], star_shaped_polygon), std::make_pair(s, x + y * image_size)));
    }
    while(!queue.empty()) {
        GrowthInfo queue_element = queue.top();
        queue.pop();
        std::pair<unsigned int, unsigned int> site_info = queue_element.second;
        unsigned int s = site_info.first;
        int x = site_info.second % image_size;
        int y = site_info.second / image_size;
        if (image_closest_site[x + y * image_size] == std::numeric_limits<unsigned int>::max()) {
            bool discard = false;
            for (int xo = -max_growth_length; (xo <= max_growth_length) && !discard; ++xo) {
                if (((x + xo) >= 0 && (x + xo) < static_cast<int>(image_size))) {
                    int xp = (x + xo);
                    for (int yo = -max_growth_length; (yo <= max_growth_length) && !discard; ++yo) {
                        if (((y + yo) >= 0 && (y + yo) < static_cast<int>(image_size))) {
                            int yp = (y + yo);
                            const Vec2d offset_length(static_cast<double>(xo), static_cast<double>(yo));
                            if (offset_length.length() <= static_cast<double>(max_growth_length)) {
                                if (image_closest_site[xp + yp * image_size] != std::numeric_limits<unsigned int>::max() && image_closest_site[xp + yp * image_size] != s) {
                                    discard = true;
                                }
                            }
                        }
                    }
                }
            }
            if (!discard) {
                image_closest_site[x + y * image_size] = s;
                for (unsigned int k = 0; k < 4; k++){
                    int xn = x + neigh_indices[k * 2];
                    if (xn >= 0 && xn < static_cast<int>(image_size)) {
                        int yn = y + neigh_indices[k * 2 + 1];
                        if (yn >= 0 && yn < static_cast<int>(image_size)) {
                            unsigned int pn = xn + yn * image_size;
                            if (image_closest_site[pn] == std::numeric_limits<unsigned int>::max()) {
                                Vec2d coordinate = get_coordinate(xn, yn, image_size);
                                queue.push(std::make_pair(get_distance(coordinate, sites[s], star_shaped_polygon), std::make_pair(s, pn)));
                            }
                        }
                    }
                }
            }
        }
    }
}

void read_points_file(std::vector<Vec2d>& points, const std::string& filename, double& min_radial_span, double& max_radial_span){
    min_radial_span = std::numeric_limits<unsigned int>::max();
    max_radial_span = -std::numeric_limits<unsigned int>::max();
    std::ifstream points_file(filename.c_str());
    if (points_file.is_open()) {
        std::string line;
        while (std::getline(points_file, line)) {
            double x_coord = std::atof(line.substr(0, line.find(",")).c_str());
            double y_coord = std::atof(line.substr(line.find(",") + 1).c_str());
            Vec2d point = Vec2d(x_coord, y_coord);
            double radial_span = point.length();
            min_radial_span = std::min(min_radial_span, radial_span);
            max_radial_span = std::max(min_radial_span, radial_span);
            points.push_back(point);
        }
        points_file.close();
    }
}

int main(int argc, char **argv) {
    assert(argc > 5);
    std::cout << "*** Parameters" << std::endl << "* distance polygon = " << argv[1] << std::endl << "* unit cell resolution = "  << argv[2] << std::endl;
    std::cout << "* output image (ppm) = " << argv[3] << std::endl << "* constrained growth length = " << argv[4] <<  std::endl << "* lattice type = " << argv[5] <<std::endl; 
    unsigned int image_size = static_cast<unsigned int>(std::atoi(argv[2]));
    assert(image_size > 0);
    int max_growth_length = atoi(argv[4]);
    assert (max_growth_length > 0);
    int lattice_type_code = atoi(argv[5]);
    LatticeType lattice_type;
    if (lattice_type_code == 0) lattice_type = DiagonalSquare;
    else if (lattice_type_code == 1) lattice_type = Triangular;
    else if (lattice_type_code == 2) lattice_type = Honeycomb;
    else return -1;
    
    std::cout << "*** Computing Voronoi growth cells..." << std::endl;
    std::vector<Vec2d> star_shaped_polygon;
    double min_radial_span, max_radial_span;
    read_points_file(star_shaped_polygon, argv[1], min_radial_span, max_radial_span);
    assert(star_shaped_polygon.size() > = 3);
    double k = 1.0; // worst case upper-bound periodic unit cell
    if (lattice_type == DiagonalSquare) k = sqrt(1.0 / 2.0);
    else if (lattice_type == Triangular) k = 1.0 / 2.0;
    else if (lattice_type == Honeycomb) k = 1.0 / 3.0;
    double periodic_offset = k * max_radial_span / min_radial_span;

    // expand to extended domain
    unsigned int image_size_ext = static_cast<unsigned int> (static_cast<double>(image_size) * (1.0 + periodic_offset * 2.0));
    std::vector<unsigned int> image_closest_site(image_size_ext * image_size_ext);
    std::fill(image_closest_site.begin(), image_closest_site.end(), std::numeric_limits<unsigned int>::max());

    // generate lattice sites
    std::vector<Vec2d> periodic_sites;
    int rings_offset = static_cast<int>(std::ceil(periodic_offset));
    int lattice_size = 2;
    if (lattice_type == Honeycomb) lattice_size = 3;
    int ini = -lattice_size * rings_offset;
    int end = lattice_size * (rings_offset + 1);
    for (int x =  ini; x < end; x++) {
        for (int y =  ini; y < end; y++) {
            Vec2d periodic_site = Vec2d(
                    (static_cast<double>(x) + 0.5) / static_cast<double>(lattice_size),
                    (static_cast<double>(y) + 0.5) / static_cast<double>(lattice_size));
            if (!(periodic_site.x >= -periodic_offset && periodic_site.x <= 1.0 + periodic_offset)) continue;
            if (!(periodic_site.y >= -periodic_offset && periodic_site.y <= 1.0 + periodic_offset)) continue;
            if (lattice_type == DiagonalSquare) {
                if (!((std::abs(x) % 2) == (std::abs(y) % 2))) continue;
            }
            if (lattice_type == Triangular || lattice_type == Honeycomb) {
                if (std::abs(y) % 2 == 1) periodic_site.x += 0.5 / static_cast<double>(lattice_size);
                periodic_site.y *= sqrt(3.0 / 4.0);
            }
            if (lattice_type == Honeycomb) {
                if (x >= 0) {
                     if (((std::abs(y) % 2 == 0) && (x % 3 == 1)) ||  ((std::abs(y) % 2 == 1) && (x % 3 == 2))) continue;
                }else{
                    if (((std::abs(y) % 2 == 1) && (std::abs(x) % 3 == 1)) || ((std::abs(y) % 2 == 0) && (std::abs(x) % 3 == 2))) continue;
                } 
            }
            periodic_sites.push_back((periodic_site + periodic_offset) / (1.0 + periodic_offset * 2.0));
        }
    }

    compute_voronoi_growth_labeling(image_size_ext, periodic_sites, image_closest_site, max_growth_length, star_shaped_polygon);

    /* export to PPM gray image (P6): a raster of Height rows, in order from top to bottom.
    Each row consists of Width pixels, in order from left to right. */
    unsigned int dim_x = image_size;
    unsigned int dim_y = image_size;
    if (lattice_type == Triangular) {
        dim_y = static_cast<unsigned int>(static_cast<double>(dim_y) * sqrt(3.0 / 4.0));
    }else if (lattice_type == Honeycomb) {
        dim_y = static_cast<unsigned int>(static_cast<double>(dim_y) * sqrt(3.0 / 4.0) * 4.0 / 3.0);
    }
    unsigned int crop_image_low = static_cast<unsigned int>(periodic_offset * static_cast<double>(image_size));
    std::ofstream ofs;
    ofs.open(argv[3], std::ios::binary);
    assert(!ofs.fail());
    ofs << "P6\n" << dim_x << " " << dim_y << "\n255\n";
    for (unsigned int y = 0; y < dim_y; y++) {
        for (unsigned int x = 0; x < dim_x ; x++) {
            unsigned char c = static_cast<unsigned char>(255);
            if (image_closest_site[(crop_image_low + x) + (crop_image_low + y) * image_size_ext] == std::numeric_limits<unsigned int>::max()){
                c = static_cast<unsigned char>(0);
            }
            ofs << c << c << c; 
        }
    }
    ofs.close();
}