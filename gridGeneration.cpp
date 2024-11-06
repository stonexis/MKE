#include <gmsh.h>
#include <vector>
#include <iostream>
#include <fstream>

int main(int argc, char **argv) {
    gmsh::initialize();
    gmsh::model::add("QuarterPlaneKirsch");
    double radius = 0.1;       // Радиус круга
    double length = 1.0;       // Размер стороны квадрата (размер четверти плоскости)
    double meshSizeNearHole = 0.01; // Размер элементов рядом с отверстием
    double meshSizeFar = 0.05;  // Размер элементов вдали от отверстия

    int rect = gmsh::model::occ::addRectangle(0, 0, 0, length, length);
    int disk = gmsh::model::occ::addDisk(0, 0, 0, radius, radius);
    gmsh::vectorpair outDimTags;
    std::vector<gmsh::vectorpair> outDimTagsMap;
   
    gmsh::model::occ::cut({{2, rect}}, {{2, disk}}, outDimTags, outDimTagsMap, -1, true, true);
    gmsh::model::occ::removeAllDuplicates();
    
    gmsh::model::occ::synchronize();
    gmsh::model::mesh::field::add("Distance", 1);
    gmsh::model::mesh::field::setNumbers(1, "NodesList", {1, 2});
    
    gmsh::model::mesh::field::add("Threshold", 2);
    gmsh::model::mesh::field::setNumber(2, "InField", 1);
    gmsh::model::mesh::field::setNumber(2, "SizeMin", meshSizeNearHole);
    gmsh::model::mesh::field::setNumber(2, "SizeMax", meshSizeFar);
    gmsh::model::mesh::field::setNumber(2, "DistMin", radius / 2);
    gmsh::model::mesh::field::setNumber(2, "DistMax", length / 2);
    gmsh::model::mesh::field::setAsBackgroundMesh(2);
    
    gmsh::option::setNumber("Mesh.RecombineAll", 0);
    gmsh::option::setNumber("Mesh.Algorithm", 6);
    gmsh::option::setNumber("Mesh.ElementOrder", 1);

    gmsh::model::mesh::generate(2);

    gmsh::write("quarter_plane_kirsch.msh");

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, paramCoords;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, paramCoords);
    std::ofstream out;
    out.open("nodes.txt");
    if (!out.is_open()) {
        throw std::runtime_error("Can't open file");
    } else {
        for (size_t i = 0; i < nodeTags.size(); ++i) {
            
            out << nodeCoords[3 * i] << ", "
            << nodeCoords[3 * i + 1] << std::endl;
    }
    out.close();
    }
    
    gmsh::fltk::run();
    gmsh::finalize();
    return 0;
}
