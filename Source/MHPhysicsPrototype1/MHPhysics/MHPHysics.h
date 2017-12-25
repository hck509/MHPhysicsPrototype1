#pragma once

#include "Math/Vector.h"

struct FMHNode
{
	FVector Position;
	float Mass;
};

struct FMHEdge
{
	int32 NodeIndices[2];
};

struct FMHTriangle
{
	int32 NodeIndices[3];
};

struct FMHMeshInfo
{
	int32 NodeIndex;
	int32 NumNodes;

	int32 EdgeIndex;
	int32 NumEdges;

	int32 TriangleIndex;
	int32 NumTriangles;

	FMHMeshInfo()
	{
		NodeIndex = 0;
		NumNodes = 0;
		EdgeIndex = 0;
		NumEdges = 0;
		TriangleIndex = 0;
		NumTriangles = 0;
	}
};

class FMHPhysics
{
public:
	FMHPhysics();
	~FMHPhysics();
	
	void GenerateFromStaticMesheActors(UWorld* World);
	FMHMeshInfo GenerateFromStaticMesh(const UStaticMesh& Mesh, const FTransform& Transform);

	void Tick(float DeltaSeconds);

	void DebugDraw(UWorld* World);

private:
	TArray<FMHNode> Nodes;
	TArray<FMHEdge> Edges;
	TArray<FMHTriangle> Triangles;
};
