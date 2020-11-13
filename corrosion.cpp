#include <dolfin.h>
#include "Poisson.h"
#include "precice/SolverInterface.hpp"

using namespace dolfin;

// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return x[0] < 0.001+DOLFIN_EPS;
  }
};

int main()
{

  auto mesh = std::make_shared<dolfin::Mesh>("Gapmesh.xml");
  
  std::vector<double> &vertices = mesh->coordinates();
  for (auto &v : vertices){
    v *= 0.001; // mesh is scaled from m to mm to better fit with the existing fracture code
  }
  
	auto V = std::make_shared<Poisson::FunctionSpace>(mesh);
	
	// Define boundary condition
	auto u0 = std::make_shared<Constant>(0.2);
	auto boundary = std::make_shared<DirichletBoundary>();
	DirichletBC bc(V, u0, boundary);

	// Define variational forms
	Poisson::BilinearForm a(V, V);
	Poisson::LinearForm L(V);
	L.f = u0;
	L.g = u0;
	
	precice::SolverInterface precice("Corrosion", "precice-config.xml", 0, 1);

  const int meshID = precice.getMeshID("Corrosion-Mesh");
  const int dataID = precice.getDataID("Gc", meshID);

  const int numberOfVertices = mesh->num_vertices();
  const int dimensions = precice.getDimensions();
  std::vector<double> writeData(numberOfVertices);
  std::vector<int> vertexIDs(numberOfVertices);
  
  precice.setMeshVertices(meshID, numberOfVertices, mesh->coordinates().data(), vertexIDs.data());

  double dt = precice.initialize();

	// solution
	Function u(V);

	// Save solution in VTK format
	File file("poisson.pvd");
	
	int i=0;

  while (precice.isCouplingOngoing()) {
		u0->operator=(0.2+i/10);
		L.f = u0;
		solve(a == L, u, bc);
		
		for (int j=0; j<numberOfVertices; j++){
		  writeData[j] = u(mesh->coordinates()[2*j], mesh->coordinates()[2*j+1]);
		}
		
    precice.writeBlockScalarData(dataID, numberOfVertices, vertexIDs.data(), writeData.data());

    dt = precice.advance(dt);
		
		file << u;
		i++;
	}
	
	precice.finalize();

	return 0;
}
