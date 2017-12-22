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

class FMHPhysics
{
public:
	FMHPhysics();
	~FMHPhysics();
	
	void GenerateFromStaticMesheActors(UWorld* World);
	void GenerateFromStaticMesh(const UStaticMesh& Mesh, const FTransform& Transform);

	void Tick(float DeltaSeconds);

	void DebugDraw(UWorld* World);

private:
	TArray<FMHNode> Nodes;
	TArray<FMHEdge> Edges;
	TArray<FMHTriangle> Triangles;
};
