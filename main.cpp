#include "include/parameterization.hpp"

int main(int argc, char** argv){
    //Mesh Setup
    string path = OBJECT_PATH + string(argv[1]);
    igl::read_triangle_mesh(path, V, F);

    //Viewer and Plugin Setup
    viewer.data().set_mesh(V, F);
    viewer.plugins.push_back(&plugin);

    //ImGuizmo Setup
    plugin.widgets.push_back(&gizmo);
    gizmo.T.block(0,3,3,1) = 0.5*(V.colwise().maxCoeff() + V.colwise().minCoeff()).transpose().cast<float>();
    const Eigen::Matrix4f T0 = gizmo.T;
    gizmo.operation = ImGuizmo::TRANSLATE;
    gizmo.callback = [&](const Eigen::Matrix4f & T){
        const Matrix4d TT = (T*T0.inverse()).cast<double>().transpose();
        viewer.data().set_vertices((V.rowwise().homogeneous()*TT).rowwise().hnormalized());
        viewer.data().compute_normals();
    };

    //ImGui Menu Setup
    plugin.widgets.push_back(&menu);

    //Workspace for Parameterization techniques
    vecXi U;
    Graph graph(V, F);
    graph.getBoundary(U);

    viewer.callback_key_pressed = [&](decltype(viewer) &,unsigned int key, int mod){
        switch(key){
            case '0':{
                viewer.data().clear();
                return true;
            }
            case '1':{
                viewer.data().set_mesh(V, F);
                return true;
            } case '2': {
                matXd Vb;
                //Slice the original mesh to obtain only the boundary vertices
                igl::slice(V, U, 1, Vb);
                viewer.data().add_points(Vb, Eigen::RowVector3d(0,1,0));
                return true;
            } case '3': {
                //Visualize the boundary edges
                matXd Vb;
                igl::slice(V, U, 1, Vb);
                matXd Vb2(Vb.rows(), 3);
                Vb2.bottomRows(Vb.rows() - 1) = Vb.topRows(Vb.rows() - 1);
                Vb2.row(0) = Vb.row(Vb.rows() - 1);
                viewer.data().add_edges(Vb, Vb2, Eigen::RowVector3d(1,0,0));
                return true;
            } case '4': {
                //Parameterize the boundary
                matXd Vb;
                igl::slice(V, U, 1, Vb);
                matXd Tb_uniform;
                auto parameters = getDiscParameters(Vb);
                getBoundaryDiscParameterization(Vb, Tb_uniform, 3 * get<3>(parameters), 0);
                matXd Tb_uniform3D;
                get3DParameterVertices(Tb_uniform, Tb_uniform3D, get<0>(parameters), get<1>(parameters), get<2>(parameters));
                viewer.data().add_points(Tb_uniform3D, Eigen::RowVector3d(0,1,1));
                return true;
            } case '5': {
                //Parameterize the boundary
                matXd Vb;
                igl::slice(V, U, 1, Vb);
                matXd Tb_chord;
                auto parameters = getDiscParameters(Vb);
                getBoundaryDiscParameterization(Vb, Tb_chord, 3 * get<3>(parameters), 1);
                matXd Tb_chord3D;
                get3DParameterVertices(Tb_chord, Tb_chord3D, get<0>(parameters), get<1>(parameters), get<2>(parameters));
                viewer.data().add_points(Tb_chord3D, Eigen::RowVector3d(1,0,0));
                return true;
            }
        }
        return false;
    };

    //Launch GUI
    viewer.launch();
}