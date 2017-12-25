#include "MHPhysics.h"
#include "DrawDebugHelpers.h"
#include "EngineUtils.h"
#include "Engine/StaticMeshActor.h"

FMHPhysics::FMHPhysics()
{

}

FMHPhysics::~FMHPhysics()
{

}

void FMHPhysics::GenerateFromStaticMesheActors(UWorld* World)
{
	for (TActorIterator<AStaticMeshActor> Iter(World); Iter; ++Iter)
	{
		AStaticMeshActor* StaticMeshActor = *Iter;

		UStaticMesh* StaticMesh = StaticMeshActor->GetStaticMeshComponent()->GetStaticMesh();
		if (StaticMesh)
		{
			GenerateFromStaticMesh(*StaticMesh, StaticMeshActor->GetActorTransform());
		}
	}
}

FMHMeshInfo FMHPhysics::GenerateFromStaticMesh(const UStaticMesh& Mesh, const FTransform& Transform)
{
	FMHMeshInfo NewMeshInfo;

	if (!ensure(Mesh.RenderData.IsValid()) || !ensure(Mesh.RenderData->LODResources.Num() > 0))
	{
		return NewMeshInfo;
	}

	NewMeshInfo.NodeIndex = Nodes.Num();
	NewMeshInfo.EdgeIndex = Edges.Num();
	NewMeshInfo.TriangleIndex = Triangles.Num();

	const int32 NodeOffset = Nodes.Num();

	const FStaticMeshLODResources& Resource = Mesh.RenderData->LODResources[0];
	
	for (uint32 VertexIndex = 0; VertexIndex < Resource.PositionVertexBuffer.GetNumVertices(); ++VertexIndex)
	{
		const FVector Position = Transform.TransformPosition(Resource.PositionVertexBuffer.VertexPosition(VertexIndex));

		Nodes.Add(FMHNode({ Position, 1.0f }));
	}

	FIndexArrayView IndexArrayView = Resource.IndexBuffer.GetArrayView();

	ensure(IndexArrayView.Num() % 3 == 0);

	for (int32 Index = 0; Index + 2 < IndexArrayView.Num(); Index += 3)
	{
		Triangles.Add(FMHTriangle({
			NodeOffset + static_cast<int32>(IndexArrayView[Index]),
			NodeOffset + static_cast<int32>(IndexArrayView[Index + 1]),
			NodeOffset + static_cast<int32>(IndexArrayView[Index + 2]) }));
	}

	NewMeshInfo.NumNodes = Nodes.Num() - NewMeshInfo.NodeIndex;
	NewMeshInfo.NumEdges = Edges.Num() - NewMeshInfo.EdgeIndex;
	NewMeshInfo.NumTriangles = Triangles.Num() - NewMeshInfo.TriangleIndex;

	return NewMeshInfo;
}

void FMHPhysics::Tick(float DeltaSeconds)
{

}

void FMHPhysics::DebugDraw(UWorld* World)
{
	for (const FMHNode& Node : Nodes)
	{
		::DrawDebugSphere(World, Node.Position, 10.0f, 3, FColor::White, false, 0);
	}

	for (const FMHTriangle& Triangle : Triangles)
	{
		::DrawDebugLine(World,
			Nodes[Triangle.NodeIndices[0]].Position,
			Nodes[Triangle.NodeIndices[1]].Position, FColor::Green, false, 0);

		::DrawDebugLine(World,
			Nodes[Triangle.NodeIndices[0]].Position,
			Nodes[Triangle.NodeIndices[1]].Position, FColor::Green, false, 0);

		::DrawDebugLine(World,
			Nodes[Triangle.NodeIndices[0]].Position,
			Nodes[Triangle.NodeIndices[2]].Position, FColor::Green, false, 0);

		::DrawDebugLine(World,
			Nodes[Triangle.NodeIndices[1]].Position,
			Nodes[Triangle.NodeIndices[2]].Position, FColor::Green, false, 0);
	}
}
