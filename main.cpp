#include "typedefs.h"
#include "lcc_comparer.h"
using namespace std;

// Functor used to display all the vertices of a given volume.
template<class LCC>
struct Display_vol_vertices : public std::unary_function<LCC, void>
{
  Display_vol_vertices(const LCC& alcc) :
    lcc(alcc),
    nb_volume(0)
  {}

  void operator() (typename LCC::Dart& d)
  {
    std::cout<<"Volume "<<++nb_volume<<" : ";
    for (typename LCC::template One_dart_per_incident_cell_range<0,3>::
           const_iterator it=lcc.template one_dart_per_incident_cell<0,3>
           (lcc.dart_handle(d)).begin(),
           itend=lcc.template one_dart_per_incident_cell<0,3>
           (lcc.dart_handle(d)).end();
         it!=itend; ++it)
    {
      std::cout << lcc.point(it) << "; ";
    }
    std::cout<<std::endl;
  }
private:
  const LCC& lcc;
  unsigned int nb_volume;
};

struct Cube
{
public:
    Dart_handle *db, *d1, *d2, *d3, *d4, *dt;
};

void display_lcc(LCC_4 lcc)
{
    std::for_each(lcc.one_dart_per_cell<3>().begin(),
                    lcc.one_dart_per_cell<3>().end(),
                    Display_vol_vertices<LCC_4>(lcc));

    lcc.display_characteristics(std::cout);
    cout<<", is_valid="<<lcc.is_valid()<<endl;
}

Point create_point(int x1, int x2, int x3, int x4)
{
    FT p[] = {x1, x2, x3, x4};
    std::vector<FT> v(p, p+4);
    return Point(4, v.begin(), v.end());
}

void print_dart_info(Dart_handle &d, LCC_4 &lcc, std::string name, ostream &stream)
{
    stream<<"Dart "<<name<<endl<<"-------"<<endl;
    stream<<"Point: "<<lcc.point(d)<<endl;
    for(unsigned int i = 0; i <= 4; i++)
    {
        Dart_handle nd = lcc.beta(d, i);
        stream<<"β"<<i<<": ";
        if (nd == lcc.null_dart_handle)
        {
            stream<<"null"<<endl;
        }
        else
        {
            stream<<lcc.point(nd)<<endl;
        }
    }
    stream<<endl;
}

void print_dart_info(Dart_handle &d, LCC_4 &lcc, std::string name)
{
    print_dart_info(d, lcc, name, cout);
}

void create_face(LCC_4 &lcc, Point &p1, Point &p2, Point &p3, Point &p4, Dart_handle *&all_dh)
{
    Dart_handle d1 = lcc.create_dart(p1);
    Dart_handle d2 = lcc.create_dart(p2);
    Dart_handle d3 = lcc.create_dart(p3);
    Dart_handle d4 = lcc.create_dart(p4);

    lcc.sew<1>(d1, d2);
    lcc.sew<1>(d2, d3);
    lcc.sew<1>(d3, d4);
    lcc.sew<1>(d4, d1);

    all_dh = new Dart_handle[4];
    all_dh[0] = d1;
    all_dh[1] = d2;
    all_dh[2] = d3;
    all_dh[3] = d4;
}

void create_cube(LCC_4 &lcc, Point &p1, Point &p2, Point &p3, Point &p4, Point &p5, Point &p6, Point &p7, Point &p8, Cube &cube)
{
    create_face(lcc, p2, p1, p4, p3, cube.db);
    create_face(lcc, p1, p2, p6, p5, cube.d1);
    create_face(lcc, p4, p1, p5, p8, cube.d2);
    create_face(lcc, p3, p4, p8, p7, cube.d3);
    create_face(lcc, p2, p3, p7, p6, cube.d4);
    create_face(lcc, p5, p6, p7, p8, cube.dt);

    // Sew bottom to side face
    lcc.sew<2>(cube.db[0], cube.d1[0]);
    lcc.sew<2>(cube.db[1], cube.d2[0]);
    lcc.sew<2>(cube.db[2], cube.d3[0]);
    lcc.sew<2>(cube.db[3], cube.d4[0]);

    // Sew top to side face
    lcc.sew<2>(cube.dt[0], cube.d1[2]);
    lcc.sew<2>(cube.dt[3], cube.d2[2]);
    lcc.sew<2>(cube.dt[2], cube.d3[2]);
    lcc.sew<2>(cube.dt[1], cube.d4[2]);

    // Sew side faces between them
    lcc.sew<2>(cube.d1[3], cube.d2[1]);
    lcc.sew<2>(cube.d2[3], cube.d3[1]);
    lcc.sew<2>(cube.d3[3], cube.d4[1]);
    lcc.sew<2>(cube.d4[3], cube.d1[1]);
}

void export_cube_to_file(Cube c, LCC_4 &lcc, string filename)
{
    ofstream outfile;
    outfile.open(filename);

    for(int i = 0; i < 4; i++)
    {
        print_dart_info(c.db[i], lcc, "db-" + std::to_string(i), outfile);
    }
    for(int i = 0; i < 4; i++)
    {
        print_dart_info(c.d1[i], lcc, "d1-" + std::to_string(i), outfile);
    }
    for(int i = 0; i < 4; i++)
    {
        print_dart_info(c.d2[i], lcc, "d2-" + std::to_string(i), outfile);
    }
    for(int i = 0; i < 4; i++)
    {
        print_dart_info(c.d3[i], lcc, "d3-" + std::to_string(i), outfile);
    }
    for(int i = 0; i < 4; i++)
    {
        print_dart_info(c.d4[i], lcc, "d4-" + std::to_string(i), outfile);
    }
    for(int i = 0; i < 4; i++)
    {
        print_dart_info(c.dt[i], lcc, "dt-" + std::to_string(i), outfile);
    }

    outfile.close();
}

std::vector<unsigned int> count_lcc_cells(LCC_4 lcc)
{
    std::vector<unsigned int> cells(6);
    for(int i = 0; i <= 5; i++)
    {
        cells[i] = i;
    }

    std::vector<unsigned int> res = lcc.count_cells(cells);

    return res;
}

int main(int argc, char *argv[])
{
    LCC_4 lcc, lcc1, lcc2;

    // Create points for the inside cube
    Point pi1 = create_point(-1,-1,-1,-1);
    Point pi2 = create_point(-1, 1,-1,-1);
    Point pi3 = create_point( 1, 1,-1,-1);
    Point pi4 = create_point( 1,-1,-1,-1);
    Point pi5 = create_point(-1,-1, 1,-1);
    Point pi6 = create_point(-1, 1, 1,-1);
    Point pi7 = create_point( 1, 1, 1,-1);
    Point pi8 = create_point( 1,-1, 1,-1);

    // Create points for the outside cube
    Point po1 = create_point(-1,-1,-1, 1);
    Point po2 = create_point(-1, 1,-1, 1);
    Point po3 = create_point( 1, 1,-1, 1);
    Point po4 = create_point( 1,-1,-1, 1);
    Point po5 = create_point(-1,-1, 1, 1);
    Point po6 = create_point(-1, 1, 1, 1);
    Point po7 = create_point( 1, 1, 1, 1);
    Point po8 = create_point( 1,-1, 1, 1);

    // Declare the tessract cubes
    Cube ci, co, c1, c2, c3, c4, cb, ct;

    // Create inside cube
    create_cube(lcc, pi1, pi2, pi3, pi4, pi5, pi6, pi7, pi8, ci);

    // Create "side" cubes
    create_cube(lcc, po1, po2, pi2, pi1, po5, po6, pi6, pi5, c1);
    create_cube(lcc, po4, po1, pi1, pi4, po8, po5, pi5, pi8, c2);
    create_cube(lcc, po3, po4, pi4, pi3, po7, po8, pi8, pi7, c3);
    create_cube(lcc, po2, po3, pi3, pi2, po6, po7, pi7, pi6, c4);

    create_cube(lcc, po1, po2, po3, po4, pi1, pi2, pi3, pi4, cb);
    create_cube(lcc, pi5, pi6, pi7, pi8, po5, po6, po7, po8, ct);

//    Dart_handle d_from_att = lcc.dart_of_attribute<0>(pi1);
//    if (d_from_att == lcc.null_dart_handle)
//    {
//        cout<<"Didn't find any darts on this point";
//    }
//    else
//    {
//        print_dart_info(d_from_att, lcc, "Found from point");
//    }

    // Create outside cube
    create_cube(lcc, po2, po1, po4, po3, po6, po5, po8, po7, co);

//    export_cube_to_file(ci, lcc, "ci-before.txt");
//    export_cube_to_file(c1, lcc, "c1-before.txt");
//    export_cube_to_file(c2, lcc, "c2-before.txt");
//    export_cube_to_file(c3, lcc, "c3-before.txt");
//    export_cube_to_file(c4, lcc, "c4-before.txt");
//    export_cube_to_file(co, lcc, "co-before.txt");
//    export_cube_to_file(cb, lcc, "cb-before.txt");
//    export_cube_to_file(ct, lcc, "ct-before.txt");

    // 3-sew the "side" and bottom-top cubes with the inside cube
    lcc.sew<3>(c1.d3[0], ci.d1[0]);
    lcc.sew<3>(c2.d3[0], ci.d2[0]);
    lcc.sew<3>(c3.d3[0], ci.d3[0]);
    lcc.sew<3>(c4.d3[0], ci.d4[0]);
    lcc.sew<3>(cb.dt[0], ci.db[0]);
    lcc.sew<3>(ct.db[0], ci.dt[0]);

    // 3-sew the "side" and bottom-top cubes with the outside cube
    lcc.sew<3>(c1.d1[0], co.d1[0]);
    lcc.sew<3>(c2.d1[0], co.d4[0]);
    lcc.sew<3>(c3.d1[0], co.d3[0]);
    lcc.sew<3>(c4.d1[0], co.d2[0]);
    lcc.sew<3>(cb.db[0], co.db[0]);
    lcc.sew<3>(ct.dt[0], co.dt[0]);

    // 3-sew the "side" cubes together
    lcc.sew<3>(c1.d2[0], c2.d4[0]);
    lcc.sew<3>(c2.d2[0], c3.d4[0]);
    lcc.sew<3>(c3.d2[0], c4.d4[0]);
    lcc.sew<3>(c4.d2[0], c1.d4[0]);

    // 3-sew the bottom cube to the "side" cubes
    lcc.sew<3>(c1.db[0], cb.d1[0]);
    lcc.sew<3>(c2.db[0], cb.d2[0]);
    lcc.sew<3>(c3.db[0], cb.d3[0]);
    lcc.sew<3>(c4.db[0], cb.d4[0]);

    // 3-sew the top cube to the "side" cubes
    lcc.sew<3>(c1.dt[0], ct.d1[2]);
    lcc.sew<3>(c2.dt[0], ct.d2[2]);
    lcc.sew<3>(c3.dt[0], ct.d3[2]);
    lcc.sew<3>(c4.dt[0], ct.d4[2]);

//    export_cube_to_file(ci, lcc, "ci-after.txt");
//    export_cube_to_file(c1, lcc, "c1-after.txt");
//    export_cube_to_file(c2, lcc, "c2-after.txt");
//    export_cube_to_file(c3, lcc, "c3-after.txt");
//    export_cube_to_file(c4, lcc, "c4-after.txt");
//    export_cube_to_file(co, lcc, "co-after.txt");
//    export_cube_to_file(cb, lcc, "cb-after.txt");
//    export_cube_to_file(ct, lcc, "ct-after.txt");

    // Create reference complex from hexahedron (for comparison)
    Dart_handle dh = lcc1.make_hexahedron(pi2, pi1, pi4, pi3, pi7, pi6, pi5, pi8);
    lcc2.make_hexahedron(pi2, pi1, pi4, pi3, pi7, pi6, pi5, pi8);

    if (lcc1.beta(dh, 1) == lcc1.beta(dh, 1))
      cout << "Are equal" << endl;

    // Print information for d1 on reference complex
    print_dart_info(dh, lcc1, "d1 on hexahedron");

    std::vector<unsigned int> res = count_lcc_cells(lcc);

    // Display information for both complexes
    cout<<"LCC0"<<endl<<"===="<<endl;
    display_lcc(lcc);

    cout<<"LCC1"<<endl<<"===="<<endl;
    display_lcc(lcc1);

    LccComparer<LCC_4> comparer;
    if (comparer.compare(lcc, lcc))
      cout << "lcc and lcc are isomorphic" << endl;
    else
      cout << "lcc and lcc are not isomorphic" << endl;

    if (comparer.compare(lcc1, lcc2))
      cout << "lcc1 and lcc2 are isomorphic" << endl;
    else
      cout << "lcc1 and lcc2 are not isomorphic" << endl;

    if (comparer.compare(lcc, lcc2))
      cout << "lcc and lcc2 are isomorphic" << endl;
    else
      cout << "lcc and lcc2 are not isomorphic" << endl;

    return EXIT_SUCCESS;
}
