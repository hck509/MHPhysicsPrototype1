#pragma once

#include "Math/Vector.h"
#include "MHPhysics.generated.h"

struct FMHNode
{
	float Mass;
	FVector Position;

	// Intermediate Data
	FVector PrevPosition;
	FVector Velocity;
	FVector Force;

	void InitNode(float InMass, const FVector& InPosition)
	{
		Mass = InMass;
		Position = InPosition;

		PrevPosition = Position;
		Velocity = FVector::ZeroVector;
		Force = FVector::ZeroVector;
	}
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

struct FMHContact
{
	int32 TriangleIndices[2];
	float Depth;
	FVector Normal;	// Triangle 2 to 1

	// Node to Plane
	struct FNodeToPlane
	{
		int32 NodeIndex;
	} NodeToPlane;
};

USTRUCT(Blueprintable)
struct FMHPhsycisSetting
{
	GENERATED_BODY()

	UPROPERTY()
	FVector Gravity;

	FMHPhsycisSetting()
	{
		Gravity.Set(0, 0, -980.0f);
	}
};

class FMHPhysics
{
public:
	FMHPhysics();
	~FMHPhysics();
	
	void GenerateFromStaticMesheActors(UWorld* World);
	FMHMeshInfo GenerateFromStaticMesh(const UStaticMesh& Mesh, const FTransform& Transform, float MeshMassInKg);

	void Tick(float DeltaSeconds);

	void DebugDraw(UWorld* World);

private:
	FMHPhsycisSetting Setting;

	TArray<FMHNode> Nodes;
	TArray<FMHEdge> Edges;
	TArray<FMHTriangle> Triangles;

	TArray<FMHContact> Contacts;
};
