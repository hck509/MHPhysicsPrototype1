#pragma once

#include "Math/Vector.h"
#include "MHPhysics.generated.h"

DECLARE_STATS_GROUP(TEXT("MHPhysics"), STATGROUP_MHP, STATCAT_Advanced);

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
struct FMHChunkDrive
{
	GENERATED_BODY()

	UPROPERTY()
	FString Name;

	UPROPERTY()
	int32 NodeIndices[2];
};

USTRUCT()
struct FMHChunkHydraulic
{
	GENERATED_BODY()

	UPROPERTY()
	FString Name;

	UPROPERTY()
	int32 NodeIndices[2];
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

	UPROPERTY()
	TArray<FMHChunkDrive> Drives;

	UPROPERTY()
	TArray<FMHChunkHydraulic> Hydraulics;

	void Clear()
	{
		Nodes.Empty();
		Edges.Empty();
		Triangles.Empty();
		Drives.Empty();
		Hydraulics.Empty();
	}

#if WITH_EDITORONLY_DATA
	bool LoadFromFbx(const FString& Filename);
#endif
};

struct FMHNode
{
	float Mass;
	FVector Position;

	bool bSelfCollision : 1;
	int32 MeshIndex;

	// Intermediate Data
	FVector PrevPosition;
	FVector Velocity;
	FVector Force;
	FBox CachedBBox;

	void InitNode(float InMass, const FVector& InPosition, bool bInSelfCollision, int32 InMeshIndex)
	{
		Mass = InMass;
		Position = InPosition;
		bSelfCollision = bInSelfCollision;
		MeshIndex = InMeshIndex;

		PrevPosition = Position;
		Velocity = FVector::ZeroVector;
		Force = FVector::ZeroVector;
		CachedBBox.IsValid = 1;
		CachedBBox.Min = CachedBBox.Max = Position;
	}

	void UpdateBBox()
	{
		CachedBBox.Min.X = FMath::Min(Position.X, PrevPosition.X);
		CachedBBox.Min.Y = FMath::Min(Position.Y, PrevPosition.Y);
		CachedBBox.Min.Z = FMath::Min(Position.Z, PrevPosition.Z);

		CachedBBox.Max.X = FMath::Max(Position.X, PrevPosition.X);
		CachedBBox.Max.Y = FMath::Max(Position.Y, PrevPosition.Y);
		CachedBBox.Max.Z = FMath::Max(Position.Z, PrevPosition.Z);
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
	int32 MeshIndex;

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

	void InitTriangle(int32 NodeIndex0, int32 NodeIndex1, int32 NodeIndex2, int32 InMeshIndex)
	{
		NodeIndices[0] = NodeIndex0;
		NodeIndices[1] = NodeIndex1;
		NodeIndices[2] = NodeIndex2;
		MeshIndex = InMeshIndex;

		CacheFlags = 0;
	}

	void UpdateBBox(const TArray<FMHNode>& Nodes);

	bool HasNodeIndex(int32 NodeIndex) const
	{
		return (NodeIndices[0] == NodeIndex) || (NodeIndices[1] == NodeIndex) || (NodeIndices[2] == NodeIndex);
	}
};

struct FMHDrive
{
	FString Name;
	int32 NodeIndices[2];
	TArray<int32> TorqueNodeIndices[2];

	float Torque;
};

struct FMHHydraulic
{
	FString Name;
	int32 NodeIndices[2];
	int32 EdgeIndex;

	float DefaultLength;
	float Length;
};

struct FMHMeshInfo
{
	int32 NodeIndex;
	int32 NumNodes;

	int32 EdgeIndex;
	int32 NumEdges;

	int32 TriangleIndex;
	int32 NumTriangles;

	int32 DriveIndex;
	int32 NumDrives;

	int32 HydraulicIndex;
	int32 NumHydraulics;

	FMHMeshInfo()
	{
		NodeIndex = 0;
		NumNodes = 0;
		EdgeIndex = 0;
		NumEdges = 0;
		TriangleIndex = 0;
		NumTriangles = 0;
		DriveIndex = 0;
		NumDrives = 0;
		HydraulicIndex = 0;
		NumHydraulics = 0;
	}
};

struct FMHMesh
{
	FMHMeshInfo Info;
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

struct FMHPhysicsProfile
{
	FMHPhysicsProfile()
		: NumStepsOnLastTick(0)
	{

	}

	void PreTick()
	{
		NumStepsOnLastTick = 0;
	}

	int32 NumStepsOnLastTick;
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
	
	void SetDriveTorque(int32 DriveIndex, float Torque);
	void SetHydraulicScale(int32 HydraulicIndex, float Scale);

	void DebugDraw(UWorld* World);

private:
	void Step(float DeltaSeconds);

	FMHPhsycisSetting Setting;

	TArray<FMHNode> Nodes;
	TArray<FMHEdge> Edges;
	TArray<FMHTriangle> Triangles;
	TArray<FMHMesh> Meshes;
	TArray<FMHDrive> Drives;
	TArray<FMHHydraulic> Hydraulics;

	// Intermediate State
	TArray<FMHContact> Contacts;

	float TickSecondLeft;

	FMHPhysicsProfile Profile;
};
