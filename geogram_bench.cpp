//#define CGAL_PROFILE
#ifndef ONLY_CGAL
#  include <geogram/basic/common.h>
#  include <geogram/basic/logger.h>
#  include <geogram/basic/command_line.h>
#  include <geogram/basic/command_line_args.h>
#  include <geogram/basic/stopwatch.h>
#  include <geogram/basic/file_system.h>
#  include <geogram/mesh/mesh.h>
#  include <geogram/mesh/mesh_io.h>
#  include <geogram/mesh/mesh_reorder.h>
#  include <geogram/delaunay/delaunay.h>
#endif
#ifndef ONLY_GEOGRAM
#  include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#  include <CGAL/IO/Triangulation_off_ostream_3.h>
#  include <CGAL/Delaunay_triangulation_3.h>
#  include <CGAL/Point_set_3.h>
#  include <CGAL/Point_set_3/IO.h>
#  include <CGAL/Timer.h>
#  include <CGAL/Random.h>
#endif
#include <algorithm>

#ifdef USE_GOOGLE_BENCHMARK
#  include <benchmark/benchmark.h>
#endif

#include <iostream>
#include <fstream>

#ifndef ONLY_CGAL
namespace {
    using namespace GEO;

   /**
    * \brief Loads points from a file.
    * \param[in] points_filename the name of the file with the points.
    *  -If the example was compiled with the Geogram library, then any
    *  mesh file handled by Geogram can be used.
    *  -if the example was compiled with Delaunay_psm (single file), then
    *  the file should be ASCII, with one point per line.
    * \param[in] dimension number of coordinates of the points.
    * \param[out] points the loaded points, in a single vector of coordinates.
    *  In the end, the number of loaded points is points.size()/dimension.
    */
    bool load_points(
        const std::string& points_filename,
        index_t dimension,
        vector<double>& points
        ) {
#ifdef GEOGRAM_PSM
        // Simple data input: one point per line, coordinates in ASCII
        LineInput input(points_filename);
        if(!input.OK()) {
            return false;
        }
        while(!input.eof() && input.get_line()) {
            input.get_fields();
            if(input.nb_fields() == dimension) {
                for(index_t c=0; c<dimension; ++c) {
                    points.push_back(input.field_as_double(c));
                }
            }
        }
#else
        // Using Geogram mesh I/O
        Mesh M;
        MeshIOFlags flags;
        flags.reset_element(MESH_FACETS);
        flags.reset_element(MESH_CELLS);
        if(!mesh_load(points_filename, M, flags)) {
            return false;
        }
        M.vertices.set_dimension(dimension);
        index_t nb_points = M.vertices.nb();
        points.resize(nb_points * dimension);
        Memory::copy(
            points.data(),
            M.vertices.point_ptr(0),
            M.vertices.nb()*dimension*sizeof(double)
        );
#endif
        return true;
    }

   /**
    * \brief Saves a Delaunay triangulation to a file.
    * \param[in] delaunay a pointer to the Delaunay triangulation.
    * \param[in] filename the name of the file to be saved.
    *  -If the example was compiled with the Geogram library, then any
    *  mesh file handled by Geogram can be used.
    *  if the example was compiled with Delaunay_psm (single file), then
    *  the points and vertices of the triangulation are output in ASCII.
    * \param[in] convex_hull_only if true, then only the triangles on the
    *  convex hull are output.
    */
    void save_Delaunay(
        Delaunay* delaunay, const std::string& filename,
        bool convex_hull_only = false
    ) {
        vector<index_t> tri2v;

        if(convex_hull_only) {

            // The convex hull can be efficiently traversed only if infinite
            // tetrahedra are kept.
            geo_assert(delaunay->keeps_infinite());

            // The convex hull can be retrieved as the finite facets
            // of the infinite cells (note: it would be also possible to
            // throw away the infinite cells and get the convex hull as
            // the facets adjacent to no cell). Here we use the infinite
            // cells to show an example with them.


            // This block is just a sanity check
            {
                for(index_t t=0; t < delaunay->nb_finite_cells(); ++t) {
                    geo_debug_assert(delaunay->cell_is_finite(t));
                }

                for(index_t t=delaunay->nb_finite_cells();
                    t < delaunay->nb_cells(); ++t) {
                    geo_debug_assert(delaunay->cell_is_infinite(t));
                }
            }

            // This iterates on the infinite cells
            for(
                index_t t = delaunay->nb_finite_cells();
                t < delaunay->nb_cells(); ++t
             ) {
                for(index_t lv=0; lv<4; ++lv) {
                    signed_index_t v = delaunay->cell_vertex(t,lv);
                    if(v != -1) {
                        tri2v.push_back(index_t(v));
                    }
                }
            }
        }

#ifdef GEOGRAM_PSM
        // Simple data output: output vertices and simplices

        Logger::out("Delaunay") << "Saving output to " << filename << std::endl;
        std::ofstream out(filename.c_str());

        out << delaunay->nb_vertices() << " vertices" << std::endl;
        for(index_t v=0; v < delaunay->nb_vertices(); ++v) {
            for(index_t c=0; c < delaunay->dimension(); ++c) {
                out << delaunay->vertex_ptr(v)[c] << " ";
            }
            out << std::endl;
        }
        if(convex_hull_only) {
            out << tri2v.size()/3 << " simplices" << std::endl;
            for(index_t t=0; t<tri2v.size()/3; ++t) {
                out << tri2v[3*t] << " "
                    << tri2v[3*t+1] << " "
                    << tri2v[3*t+2] << std::endl;
            }
        } else {
            out << delaunay->nb_cells() << " simplices" << std::endl;
            for(index_t t=0; t<delaunay->nb_cells(); ++t) {
                for(index_t lv=0; lv<delaunay->cell_size(); ++lv) {
                    out << delaunay->cell_vertex(t,lv) << " ";
                }
                out << std::endl;
            }
        }

#else
        // Using Geogram mesh I/O: copy Delaunay into a Geogram
        // mesh and save it to disk.

        Mesh M_out;
        vector<double> pts(delaunay->nb_vertices() * 3);
        for(index_t v = 0; v < delaunay->nb_vertices(); ++v) {
            pts[3 * v] = delaunay->vertex_ptr(v)[0];
            pts[3 * v + 1] = delaunay->vertex_ptr(v)[1];
            pts[3 * v + 2] =
                (delaunay->dimension() >= 3) ? delaunay->vertex_ptr(v)[2] : 0.0;
        }

        if(convex_hull_only) {
            M_out.facets.assign_triangle_mesh(3, pts, tri2v, true);
        } else if(delaunay->dimension() == 3) {
            vector<index_t> tet2v(delaunay->nb_cells() * 4);
            for(index_t t = 0; t < delaunay->nb_cells(); ++t) {
            tet2v[4 * t] = index_t(delaunay->cell_vertex(t, 0));
            tet2v[4 * t + 1] = index_t(delaunay->cell_vertex(t, 1));
            tet2v[4 * t + 2] = index_t(delaunay->cell_vertex(t, 2));
            tet2v[4 * t + 3] = index_t(delaunay->cell_vertex(t, 3));
            }
            M_out.cells.assign_tet_mesh(3, pts, tet2v, true);
        } else if(delaunay->dimension() == 2) {
            tri2v.resize(delaunay->nb_cells() * 3);
            for(index_t t = 0; t < delaunay->nb_cells(); ++t) {
                tri2v[3 * t] = index_t(delaunay->cell_vertex(t, 0));
                tri2v[3 * t + 1] = index_t(delaunay->cell_vertex(t, 1));
                tri2v[3 * t + 2] = index_t(delaunay->cell_vertex(t, 2));
            }
            M_out.facets.assign_triangle_mesh(3, pts, tri2v, true);
        }
        M_out.show_stats();

        Logger::div("Saving the result");
        MeshIOFlags flags;
        flags.set_element(MESH_FACETS);
        flags.set_element(MESH_CELLS);
        mesh_save(M_out, filename, flags);
#endif
    }

}
#endif // !ONLY_CGAL

#ifndef ONLY_GEOGRAM
CGAL::Timer timer;
#endif // ! ONLY_GEOGRAM

// global variables used by bench_geogram and bench_cgal
int argc;
char** argv;


#ifndef ONLY_CGAL
static void bench_geogram(
#if USE_GOOGLE_BENCHMARK
                       benchmark::State& state
#else
                       std::array<int, 1> state = { 0 }
#endif
)
{
    using namespace GEO;
    // Needs to be called once.
    GEO::initialize();

    {
        std::vector<std::string> filenames;

        CmdLine::import_arg_group("standard");
        CmdLine::import_arg_group("algo");
        CmdLine::set_arg("algo:delaunay","default");
        CmdLine::set_arg("sys:multithread","false");

        if(
            !CmdLine::parse(
                argc, argv, filenames, "pointsfile <outputfile|none>"
            )
        ) {
            ::exit(1);
        }


        std::string points_filename = filenames[0];

        std::string output_filename =
            filenames.size() >= 2 ? filenames[1] : std::string("none");

        bool output = (output_filename != "none");

        Logger::div("Data I/O");

        Logger::out("I/O") << "Output = " << output_filename << std::endl;

        vector<double> points;

        if(!load_points(points_filename, 3, points)) {
            Logger::err("Delaunay") << "Could not load points" << std::endl;
            ::exit(1);
        }
        index_t nb_points = points.size() / 3;

        Logger::out("Delaunay")
            << "Loaded " << nb_points << " points" << std::endl;

        double time = 0.0;
        for([[maybe_unused]] auto _ : state) {
#ifndef ONLY_GEOGRAM
          timer.start();                                                              // <---- Start timer
#endif // ! ONLY_GEOGRAM
          Stopwatch Wdel("Delaunay");
          
          // Note: To create a parallel Delaunay 3D, one can use directly:
          // Delaunay_var delaunay = Delaunay::create(3,"PDEL") instead
          // of the line below (that uses the command line to select the
          // implementation of Delaunay).
          Delaunay_var delaunay = Delaunay::create(3, "BDEL");

          // Note: this does not transfer ownership of memory, caller
          // is still responsible of the memory of the points (here the
          // vector<double>). No memory is copied, Delaunay just keeps
          // a pointer.
          delaunay->set_vertices(nb_points, points.data());

          time = Wdel.elapsed_time();
#ifndef ONLY_GEOGRAM
          timer.stop();                                                           //<------ STOP timer
#endif // ! ONLY_GEOGRAM

          Logger::out("Delaunay") << delaunay->nb_cells() << " tetrahedra"
                                  << std::endl;

          Logger::out("Delaunay") << double(delaunay->nb_cells()) / time
                                  << " tetrahedra / second"
                                  << std::endl;
          if(output)
            save_Delaunay(delaunay, output_filename);
        }
    }
#ifndef ONLY_GEOGRAM
  std::cout<<"CGAL Timer says : "<<timer.time()<<"s."<<std::endl;
#endif // ! ONLY_GEOGRAM
}
#  ifdef USE_GOOGLE_BENCHMARK
     BENCHMARK(bench_geogram)->Unit(benchmark::kMillisecond);
#endif
#endif // ONLY_CGAL

#ifndef ONLY_GEOGRAM
static void bench_cgal(
#if USE_GOOGLE_BENCHMARK
                       benchmark::State& state
#else
                       std::array<int, 1> state = { 0 }
#endif
)
{
    std::string points_filename = argv[1];
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Delaunay_triangulation_3<K>      Triangulation;
    typedef Triangulation::Point          Point;
    typedef CGAL::Point_set_3<Point> Point_set;

    CGAL::get_default_random() = CGAL::Random(0);

    Point_set ps;
    std::ifstream in (points_filename);
    in >> ps;
    in.close();
    std::cout<<ps.size()<<" points were loaded."<<std::endl;
    for([[maybe_unused]] auto _ : state)
    {
      timer.reset();
      timer.start();
      Triangulation tr ( ps.points().begin(), ps.points().end());
      CGAL_USE(tr);
      timer.stop();
    }
    std::cout<<" CGAL Timer says :"<<timer.time()<<"s."<<std::endl;
    // std::ofstream off_stream("cgal_out.off");
    // CGAL::export_triangulation_3_to_off(off_stream, tr);
}
#  ifdef USE_GOOGLE_BENCHMARK
     BENCHMARK(bench_cgal)->Unit(benchmark::kMillisecond);
#endif
#endif // !ONLY_GEOGRAM

int main(int argc, char** argv) {
#ifdef USE_GOOGLE_BENCHMARK
  benchmark::Initialize(&argc, argv);
#endif
  ::argc = argc;
  ::argv = argv;
#ifndef ONLY_CGAL
  /*--------------
  | GEOGRAM PART |
  ---------------*/
#  ifndef USE_GOOGLE_BENCHMARK
  bench_geogram();
#  endif
#endif // ONLY_CGAL

#ifndef ONLY_GEOGRAM
  CGAL_USE(argc);
  CGAL_USE(argv);
  /*----------
  | CGAL PART |
  -------------*/
#  ifndef USE_GOOGLE_BENCHMARK
  bench_cgal();
#  endif
#endif // !ONLY_GEOGRAM
#if USE_GOOGLE_BENCHMARK
  benchmark::RunSpecifiedBenchmarks();
#endif
  return 0;
}
