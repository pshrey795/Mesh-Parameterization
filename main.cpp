#include "include/graph.hpp"

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
    matXd Vb;
    // igl::slice(V, U, 1, Vb);
    Vb.resize(U.rows(), 3);
    for(int i = 0; i < U.rows(); i++){
        Vb.row(i) = V.row(U(i));
    }

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
                viewer.data().add_points(Vb, Eigen::RowVector3d(0,1,0));
                return true;
            } case '3': {
                //Visualize the boundary edges
                matXd Vb2(Vb.rows(), 3);
                Vb2.bottomRows(Vb.rows() - 1) = Vb.topRows(Vb.rows() - 1);
                Vb2.row(0) = Vb.row(Vb.rows() - 1);
                viewer.data().add_edges(Vb, Vb2, Eigen::RowVector3d(1,0,0));
                return true;
            } case '4': {
                //Parameterize the boundary
                matXd Tb_uniform;
                auto parameters = getDiscParameters(Vb);
                getBoundaryDiscParameterization(Vb, Tb_uniform, 3 * get<3>(parameters), 0);
                matXd Tb_uniform3D;
                get3DParameterVertices(Tb_uniform, Tb_uniform3D, get<0>(parameters), get<1>(parameters), get<2>(parameters));
                viewer.data().add_points(Tb_uniform3D, Eigen::RowVector3d(0,1,1));
                viewer.data().add_edges(Vb, Tb_uniform3D, Eigen::RowVector3d(1,0,0));
                return true;
            } case '5': {
                //Parameterize the boundary
                matXd Tb_chord;
                auto parameters = getDiscParameters(Vb);
                getBoundaryDiscParameterization(Vb, Tb_chord, 3 * get<3>(parameters), 1);
                matXd Tb_chord3D;
                get3DParameterVertices(Tb_chord, Tb_chord3D, get<0>(parameters), get<1>(parameters), get<2>(parameters));
                viewer.data().add_points(Tb_chord3D, Eigen::RowVector3d(1,0,0));
                viewer.data().add_edges(Vb, Tb_chord3D, Eigen::RowVector3d(1,0,0));
                return true;
            } case '6': {
                //Obtain parameterization of internal vertices 
                matXd Tb_uniform;
                auto parameters = getDiscParameters(Vb);
                getBoundaryDiscParameterization(Vb, Tb_uniform, 3 * get<3>(parameters), 0);
                matXd T_uniform;
                graph.getInternalParameterization(Tb_uniform, T_uniform, 0);
                matXd T_uniform3D;
                get3DParameterVertices(T_uniform, T_uniform3D, get<0>(parameters), get<1>(parameters), get<2>(parameters));
                matXd V_r, T_r;
                V_r.resize(10,3);
                T_r.resize(10,3);
                for(int i = 0; i < 10; i++){
                    int k = V.rows() * (rand() / (RAND_MAX + 1.0));
                    V_r.row(i) = V.row(k);
                    T_r.row(i) = T_uniform3D.row(k);
                }
                viewer.data().add_edges(V_r, T_r, Eigen::RowVector3d(1,0,0));
                // viewer.data().set_mesh(T_uniform3D, F);
                // matXd T_uniform_b;
                // igl::slice(T_uniform3D, U, 1, T_uniform_b);
                // viewer.data().add_points(T_uniform_b, Eigen::RowVector3d(1,0,0));
                return true;
            } case '7': {
                //Obtain parameterization of internal vertices
                matXd Tb_chord;
                auto parameters = getDiscParameters(Vb);
                getBoundaryDiscParameterization(Vb, Tb_chord, 3 * get<3>(parameters), 1);
                matXd T_chord;
                graph.getInternalParameterization(Tb_chord, T_chord, 0);
                matXd T_chord3D;
                get3DParameterVertices(T_chord, T_chord3D, get<0>(parameters), get<1>(parameters), get<2>(parameters));
                viewer.data().set_mesh(T_chord3D, F);
                // matXd T_uniform_b;
                // igl::slice(T_chord3D, U, 1, T_uniform_b);
                // viewer.data().add_points(T_uniform_b, Eigen::RowVector3d(1,0,0));
                return true;
            }
        }
        return false;
    };

    //Launch GUI
    viewer.launch();
}