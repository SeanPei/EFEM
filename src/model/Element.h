#ifndef EFEM_MODEL_ELEMENT_H
#define EFEM_MODEL_ELEMENT_H

#include "core/Object.h"
#include "io/Logger.h"
#include "math/Math.h"
#include "model/Node.h"
#include <vector>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkIdList.h>
#include <vtkMassProperties.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkUnstructuredGrid.h>

namespace EFEM {
class Element : public Object {
public:
    std::vector<int> node_indexs;
    std::vector<std::shared_ptr<Node>> nodes;
    int material_id;
    std::unordered_map<std::string, std::shared_ptr<Object>> extra_info;

    void init()
    {
        m_vtk_element_data = vtkSmartPointer<vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkPoints> vtk_element_points = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkIdList> vtk_cell_point_lists = vtkSmartPointer<vtkIdList>::New();
        for (int i = 0; i < nodes.size(); i++) {
            vtk_cell_point_lists->InsertNextId(i);
            Eigen::Vector3d position = nodes[i]->position + nodes[i]->dx;
            vtk_element_points->InsertNextPoint(position.data());
        }
        VTKCellType cell_type;
        switch (nodes.size()) {
        case 4:
            cell_type = VTKCellType::VTK_TETRA;
            break;
        case 8:
            cell_type = VTKCellType::VTK_HEXAHEDRON;
            break;
        case 10:
            cell_type = VTKCellType::VTK_QUADRATIC_TETRA;
            break;
        case 20:
            cell_type = VTKCellType::VTK_QUADRATIC_HEXAHEDRON;
            break;
        }
        m_vtk_element_data->SetPoints(vtk_element_points);
        m_vtk_element_data->InsertNextCell(cell_type, vtk_cell_point_lists);

        vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
        surfaceFilter->SetInputData(m_vtk_element_data);
        surfaceFilter->Update();
        m_vtk_poly_data = surfaceFilter->GetOutput();

        vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
        triangleFilter->SetInputData(m_vtk_poly_data);
        triangleFilter->Update();
        m_vtk_poly_data = triangleFilter->GetOutput();

        m_surface_nodes.clear();
        int surf_point_num = m_vtk_poly_data->GetNumberOfPoints();
        for (int i = 0; i < surf_point_num; i++) {
            double* poly_point_coord = m_vtk_poly_data->GetPoint(i);
            int poly_point_index = m_vtk_element_data->FindPoint(poly_point_coord);
            m_surface_nodes.push_back(nodes[poly_point_index]);
        }
    }

    double volume()
    {
        if (m_surface_nodes.size() <= 0) {
            Logger::info("element not initial!");
            init();
        }
        vtkSmartPointer<vtkPoints> vtk_surface_points = vtkSmartPointer<vtkPoints>::New();
        for (int i = 0; i < m_surface_nodes.size(); i++) {
            Eigen::Vector3d position = m_surface_nodes[i]->position + m_surface_nodes[i]->dx;
            vtk_surface_points->InsertNextPoint(position.data());
        }
        m_vtk_poly_data->SetPoints(vtk_surface_points);

        vtkSmartPointer<vtkMassProperties> poly_properties = vtkSmartPointer<vtkMassProperties>::New();
        poly_properties->SetInputData(m_vtk_poly_data);
        poly_properties->Update();
        return poly_properties->GetVolumeProjected();
    }

private:
    vtkSmartPointer<vtkUnstructuredGrid> m_vtk_element_data;
    vtkSmartPointer<vtkPolyData> m_vtk_poly_data;
    std::vector<std::shared_ptr<Node>> m_surface_nodes;
};
}

#endif