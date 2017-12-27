#pragma once

#include "Math/Vector.h"
#include "MHPhysics.generated.h"

DECLARE_STATS_GROUP(TEXT("MHPhysics"), STATGROUP_MHPhysics, STATCAT_Advanced);

USTRUCT()
struct FMHChunkNode
{
	GENERATED_BODY()

	UPROPERTY()
	FVector Position;
};

USTRUCT()
struct FMHChunkEdge
{
	GENERATED_BODY()

	UPROPERTY()
	int32 NodeIndices[2];
};

USTRUCT()
struct FMHChunkTriangle
{
	GENERATED_BODY()

	UPROPERTY()
	int32 NodeIndices[3];
};

USTRUCT()
struct FMHChunk
{
	GENERATED_BODY()

	UPROPERTY()
	TArray<FMHChunkNode> Nodes;

	UPROPERTY()
	TArray<FMHChunkEdge> Edges;

	UPROPERTY()
	TArray<FMHChunkTriangle> Triangles;

	void Clear()
	{
		Nodes.Empty();
		Edges.Empty();
		Triangles.Empty();
	}

#if WITH_EDITORONLY_DATA
	bool LoadFromFbx(const FString& Filename);
#endif
};

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
	float SpringK;
	float SpringD;

	float DefaultLength;
};

struct FMHTriangle
{
	int32 NodeIndices[3];

	enum ECacheBitFlags
	{
		ValidPlane = (1 << 0),
		ValidPrevPlane = (1 << 1)
	};

	// Cache
	mutable uint32 CacheFlags;
	mutable FPlane CachedPlane;
	mutable FPlane CachedPrevPlane;
	mutable FBox CachedBBox;	// Axis-aligned bounding box from prev position to current position

	FMHTriangle()
	{
		CacheFlags = 0;
	}

	void UpdateBBox(const TArray<FMHNode>& Nodes)
	{
		CachedBBox.Min.X = FMath::Min3(Nodes[NodeIndices[0]].Position.X, Nodes[NodeIndices[1]].Position.X, Nodes[NodeIndices[2]].Position.X);
		CachedBBox.Min.Y = FMath::Min3(Nodes[NodeIndices[0]].Position.Y, Nodes[NodeIndices[1]].Position.Y, Nodes[NodeIndices[2]].Position.Y);
		CachedBBox.Min.Z = FMath::Min3(Nodes[NodeIndices[0]].Position.Z, Nodes[NodeIndices[1]].Position.Z, Nodes[NodeIndices[2]].Position.Z);

		CachedBBox.Max.X = FMath::Max3(Nodes[NodeIndices[0]].Position.X, Nodes[NodeIndices[1]].Position.X, Nodes[NodeIndices[2]].Position.X);
		CachedBBox.Max.Y = FMath::Max3(Nodes[NodeIndices[0]].Position.Y, Nodes[NodeIndices[1]].Position.Y, Nodes[NodeIndices[2]].Position.Y);
		CachedBBox.Max.Z = FMath::Max3(Nodes[NodeIndices[0]].Position.Z, Nodes[NodeIndices[1]].Position.Z, Nodes[NodeIndices[2]].Position.Z);

		CachedBBox.Min.X = FMath::Min3(FMath::Min(CachedBBox.Min.X, Nodes[NodeIndices[0]].PrevPosition.X), Nodes[NodeIndices[1]].PrevPosition.X, Nodes[NodeIndices[2]].PrevPosition.X);
		CachedBBox.Min.Y = FMath::Min3(FMath::Min(CachedBBox.Min.Y, Nodes[NodeIndices[0]].PrevPosition.Y), Nodes[NodeIndices[1]].PrevPosition.Y, Nodes[NodeIndices[2]].PrevPosition.Y);
		CachedBBox.Min.Z = FMath::Min3(FMath::Min(CachedBBox.Min.Z, Nodes[NodeIndices[0]].PrevPosition.Z), Nodes[NodeIndices[1]].PrevPosition.Z, Nodes[NodeIndices[2]].PrevPosition.Z);

		CachedBBox.Max.X = FMath::Max3(FMath::Max(CachedBBox.Max.X, Nodes[NodeIndices[0]].PrevPosition.X), Nodes[NodeIndices[1]].PrevPosition.X, Nodes[NodeIndices[2]].PrevPosition.X);
		CachedBBox.Max.Y = FMath::Max3(FMath::Max(CachedBBox.Max.Y, Nodes[NodeIndices[0]].PrevPosition.Y), Nodes[NodeIndices[1]].PrevPosition.Y, Nodes[NodeIndices[2]].PrevPosition.Y);
		CachedBBox.Max.Z = FMath::Max3(FMath::Max(CachedBBox.Max.Z, Nodes[NodeIndices[0]].PrevPosition.Z), Nodes[NodeIndices[1]].PrevPosition.Z, Nodes[NodeIndices[2]].PrevPosition.Z);
	}

	bool HasNodeIndex(int32 NodeIndex) const
	{
		return (NodeIndices[0] == NodeIndex) || (NodeIndices[1] == NodeIndex) || (NodeIndices[2] == NodeIndex);
	}
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

	enum class EType
	{
		NodeToTriangle
	} Type;

	// Node to Triangle. TriangleIndices[1] is triangle side.
	struct FNodeToTriangle
	{
		int32 NodeIndex;
	} NodeToTriangle;
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
	FMHMeshInfo GenerateFromStaticMesh(const UStaticMesh& Mesh, const FTransform& Transform, float MeshMassInKg, float SpringK, float SpringD);
	FMHMeshInfo GenerateFromChunk(const FMHChunk& Chunk, const FTransform& Transform, float MeshMassInKg, float SpringK, float SpringD);

	void Tick(float DeltaSeconds);

	const FMHNode* FindNode(int32 NodeIndex) const;
	const FMHTriangle* FindTriangle(int32 TriangleIndex) const;

	void DebugDraw(UWorld* World);

private:
	FMHPhsycisSetting Setting;

	TArray<FMHNode> Nodes;
	TArray<FMHEdge> Edges;
	TArray<FMHTriangle> Triangles;

	TArray<FMHContact> Contacts;
};
