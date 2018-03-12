#include "MHPhysics.h"
#include "DrawDebugHelpers.h"
#include "EngineUtils.h"
#include "Engine/StaticMeshActor.h"
#include "Engine.h"

#if WITH_EDITOR
#include "FbxImporter.h"
#include "Misc/Paths.h"
#endif

static TAutoConsoleVariable<float> CVarMHPhysicsFPS(
	TEXT("mhp.fps"),
	2000.0f,
	TEXT("")
);

static TAutoConsoleVariable<float> CVarMHPhysicsSpeed(
	TEXT("mhp.speed"),
	1.0f,
	TEXT("Speed multiplier")
);

static TAutoConsoleVariable<int32> CVarMHPhysicsDraw(
	TEXT("mhp.draw"),
	1,
	TEXT("")
);

static TAutoConsoleVariable<int32> CVarMHPhysicsEnableStaticNodeCollision(
	TEXT("mhp.enableStaticNodeCollision"),
	0,
	TEXT("")
);

static TAutoConsoleVariable<float> CVarMHPhysicsFrictionCoeff(
	TEXT("mhp.frictionCoeff"),
	0.5f,
	TEXT("")
);

static TAutoConsoleVariable<float> CVarMHPMaxSpeed(
	TEXT("mhp.maxSpeed"),
	1000.0f,
	TEXT("cm/s")
);

DECLARE_CYCLE_STAT(TEXT("MHP Tick"), STAT_MHPhysicsTick, STATGROUP_MHP);
DECLARE_CYCLE_STAT(TEXT("MHP Step"), STAT_MHPhysicsStep, STATGROUP_MHP);
DECLARE_CYCLE_STAT(TEXT("MHP CollisionDetect"), STAT_CollisionDetect, STATGROUP_MHP);
DECLARE_CYCLE_STAT(TEXT("MHP CD Triangle to Triangle"), STAT_CollisionDetectTriangleToTriangle, STATGROUP_MHP);
DECLARE_CYCLE_STAT(TEXT("MHP CD Node to Triangle"), STAT_CollisionDetectNodeToTriangle, STATGROUP_MHP);
DECLARE_CYCLE_STAT(TEXT("MHP Gravity"), STAT_Gravity, STATGROUP_MHP);
DECLARE_CYCLE_STAT(TEXT("MHP Edge Spring"), STAT_EdgeSpring, STATGROUP_MHP);
DECLARE_CYCLE_STAT(TEXT("MHP CollisionResolve"), STAT_CollisionResolve, STATGROUP_MHP);
DECLARE_CYCLE_STAT(TEXT("MHP Integration"), STAT_Integration, STATGROUP_MHP);
DECLARE_CYCLE_STAT(TEXT("MHP UpdateNodeBBox"), STAT_UpdateNodeBBox, STATGROUP_MHP);
DECLARE_CYCLE_STAT(TEXT("MHP UpdateTriangleBBox"), STAT_UpdateTriangleBBox, STATGROUP_MHP);
DECLARE_CYCLE_STAT(TEXT("MHP Draw"), STAT_Draw, STATGROUP_MHP);


// Intersection Code from internet :)
static bool _RayTriangleIntersect(
	const FVector &Origin, const FVector &Dir,
	const FVector &V0, const FVector &V1, const FVector &V2, const FVector& N,
	float &OutDepth)
{
	// check if ray and plane are parallel?
	float NoR = N | Dir;

	if (FMath::Abs(NoR) < SMALL_NUMBER)
	{
		// they are parallel so they don't intersect!
		return false;
	}

	// compute d parameter
	float d = N | V0;

	// compute t
	OutDepth = (d - (N | Origin)) / NoR;
	
	// check if the triangle is in behind the ray
	if (OutDepth < 0)
	{
		// the triangle is behind 
		return false;
	}

	// compute the intersection point
	FVector P = Origin + OutDepth * Dir;

	// Step 2: inside-outside test
	FVector C; // vector perpendicular to triangle's plane 

	// edge 0
	FVector edge0 = V1 - V0;
	FVector vp0 = P - V0;
	C = edge0 ^ vp0;
	if ((N | C) < 0)
	{
		// P is on the right side 
		return false;
	}

	// edge 1
	FVector edge1 = V2 - V1;
	FVector vp1 = P - V1;
	C = edge1 ^ vp1;
	if ((N | C) < 0)
	{
		// P is on the right side 
		return false;
	}

	// edge 2
	FVector edge2 = V0 - V2;
	FVector vp2 = P - V2;
	C = edge2 ^ vp2;
	if ((N | C) < 0)
	{
		// P is on the right side; 
		return false;
	}

	// this ray hits the triangle
	return true; 
}


FMHPhysics::FMHPhysics()
{
	TickSecondLeft = 0;
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
			GenerateFromStaticMesh(*StaticMesh, StaticMeshActor->GetActorTransform(), 0.0f, 0.0f, 0.0f);
		}
	}
}

FMHMeshInfo FMHPhysics::GenerateFromStaticMesh(const UStaticMesh& Mesh, const FTransform& Transform, float MeshMassInKg, float SpringK, float SpringD)
{
	FMHMeshInfo NewMeshInfo;

	if (!ensure(Mesh.RenderData.IsValid()) || !ensure(Mesh.RenderData->LODResources.Num() > 0))
	{
		return NewMeshInfo;
	}

	const int32 MeshIndex = Meshes.Num();

	NewMeshInfo.NodeIndex = Nodes.Num();
	NewMeshInfo.EdgeIndex = Edges.Num();
	NewMeshInfo.TriangleIndex = Triangles.Num();

	const int32 NodeOffset = Nodes.Num();

	const FStaticMeshLODResources& Resource = Mesh.RenderData->LODResources[0];
	const uint32 NumVertices = Resource.PositionVertexBuffer.GetNumVertices();
	const float MassPerNode = NumVertices > 0 ? MeshMassInKg / NumVertices : 0;
	
	for (uint32 VertexIndex = 0; VertexIndex < NumVertices; ++VertexIndex)
	{
		const FVector Position = Transform.TransformPosition(Resource.PositionVertexBuffer.VertexPosition(VertexIndex));

		FMHNode NewNode;
		NewNode.InitNode(MassPerNode, Position, false, MeshIndex);

		Nodes.Add(NewNode);
	}

	FIndexArrayView IndexArrayView = Resource.IndexBuffer.GetArrayView();

	ensure(IndexArrayView.Num() % 3 == 0);

	Triangles.Reserve(Triangles.Num() + IndexArrayView.Num() / 3);

	for (int32 Index = 0; Index + 2 < IndexArrayView.Num(); Index += 3)
	{
		Triangles.AddDefaulted(1);

		FMHTriangle& NewTriangle = Triangles.Last();

		NewTriangle.InitTriangle(
			NodeOffset + static_cast<int32>(IndexArrayView[Index]),
			NodeOffset + static_cast<int32>(IndexArrayView[Index + 2]),
			NodeOffset + static_cast<int32>(IndexArrayView[Index + 1]),
			MeshIndex);

		NewTriangle.UpdateBBox(Nodes);

		float Distances[3] = {
			(Nodes[NewTriangle.NodeIndices[0]].Position - Nodes[NewTriangle.NodeIndices[1]].Position).Size(),
			(Nodes[NewTriangle.NodeIndices[0]].Position - Nodes[NewTriangle.NodeIndices[2]].Position).Size(),
			(Nodes[NewTriangle.NodeIndices[1]].Position - Nodes[NewTriangle.NodeIndices[2]].Position).Size()
		};

		Edges.Add(FMHEdge({ NewTriangle.NodeIndices[0], NewTriangle.NodeIndices[1], SpringK, SpringD, Distances[0] }));
		Edges.Add(FMHEdge({ NewTriangle.NodeIndices[0], NewTriangle.NodeIndices[2], SpringK, SpringD, Distances[1] }));
		Edges.Add(FMHEdge({ NewTriangle.NodeIndices[1], NewTriangle.NodeIndices[2], SpringK, SpringD, Distances[2] }));
	}

	NewMeshInfo.NumNodes = Nodes.Num() - NewMeshInfo.NodeIndex;
	NewMeshInfo.NumEdges = Edges.Num() - NewMeshInfo.EdgeIndex;
	NewMeshInfo.NumTriangles = Triangles.Num() - NewMeshInfo.TriangleIndex;

	FMHMesh NewMesh;
	NewMesh.Info = NewMeshInfo;

	Meshes.Add(NewMesh);

	return NewMeshInfo;
}

FMHMeshInfo FMHPhysics::GenerateFromChunk(const FMHChunk& Chunk, const FTransform& Transform, float MeshMassInKg, float SpringK, float SpringD)
{
	FMHMeshInfo NewMeshInfo;

	const int32 MeshIndex = Meshes.Num();

	NewMeshInfo.NodeIndex = Nodes.Num();
	NewMeshInfo.EdgeIndex = Edges.Num();
	NewMeshInfo.TriangleIndex = Triangles.Num();
	NewMeshInfo.DriveIndex = Drives.Num();
	NewMeshInfo.HydraulicIndex = Hydraulics.Num();

	const int32 NodeOffset = Nodes.Num();
	const int32 NumNodes = Chunk.Nodes.Num();
	const float MassPerNode = NumNodes > 0 ? MeshMassInKg / NumNodes : 0;

	for (int32 i = 0; i < NumNodes; ++i)
	{
		FMHNode NewNode;
		NewNode.InitNode(MassPerNode, Transform.TransformPosition(Chunk.Nodes[i].Position), false, MeshIndex);

		Nodes.Add(NewNode);
	}

	const int32 NumEdges = Chunk.Edges.Num();
	for (int32 i = 0; i < NumEdges; ++i)
	{
		const FMHChunkEdge& Edge = Chunk.Edges[i];
		int32 NodeIndices[2] = {
			NodeOffset + Edge.NodeIndices[0],
			NodeOffset + Edge.NodeIndices[1]
		};

		const float Distance = (Nodes[NodeIndices[0]].Position - Nodes[NodeIndices[1]].Position).Size();

		Edges.Add(FMHEdge({ NodeIndices[0], NodeIndices[1], SpringK, SpringD, Distance }));
	}

	const int32 NumTriangles = Chunk.Triangles.Num();
	for (int32 i = 0; i < NumTriangles; ++i)
	{
		Triangles.AddDefaulted(1);

		FMHTriangle& NewTriangle = Triangles.Last();

		NewTriangle.InitTriangle(
			NodeOffset + Chunk.Triangles[i].NodeIndices[0],
			NodeOffset + Chunk.Triangles[i].NodeIndices[1],
			NodeOffset + Chunk.Triangles[i].NodeIndices[2],
			MeshIndex);

		NewTriangle.UpdateBBox(Nodes);
	}

	for (const FMHChunkDrive& Drive : Chunk.Drives)
	{
		Drives.AddDefaulted(1);

		FMHDrive& NewDrive = Drives.Last();
		NewDrive.Name = Drive.Name;
		NewDrive.NodeIndices[0] = NodeOffset + Drive.NodeIndices[0];
		NewDrive.NodeIndices[1] = NodeOffset + Drive.NodeIndices[1];
		NewDrive.Torque = 0;

		TArray<int32> Drive0EdgedNodes;
		TArray<int32> Drive1EdgedNodes;

		// Find neighboring edges
		const int32 NumEdges = Chunk.Edges.Num();
		for (int32 i = 0; i < NumEdges; ++i)
		{
			const FMHEdge& Edge = Edges[NewMeshInfo.EdgeIndex + i];

			if (Edge.NodeIndices[0] == NewDrive.NodeIndices[0])
			{
				Drive0EdgedNodes.AddUnique(Edge.NodeIndices[1]);
			}

			if (Edge.NodeIndices[1] == NewDrive.NodeIndices[0])
			{
				Drive0EdgedNodes.AddUnique(Edge.NodeIndices[0]);
			}

			if (Edge.NodeIndices[0] == NewDrive.NodeIndices[1])
			{
				Drive1EdgedNodes.AddUnique(Edge.NodeIndices[1]);
			}

			if (Edge.NodeIndices[1] == NewDrive.NodeIndices[1])
			{
				Drive1EdgedNodes.AddUnique(Edge.NodeIndices[0]);
			}
		}

		const FVector NodePositions[2] = {
			Nodes[NewDrive.NodeIndices[0]].Position,
			Nodes[NewDrive.NodeIndices[1]].Position
		};

		const FVector DriveAxis = NodePositions[1] - NodePositions[0];

		for (int32 i : Drive0EdgedNodes)
		{
			if (Drive1EdgedNodes.Contains(i))
			{
				// See if which side?
				const float DistanceTo0 = FMath::Abs((Nodes[i].Position - NodePositions[0]) | DriveAxis);
				const float DistanceTo1 = FMath::Abs((Nodes[i].Position - NodePositions[1]) | DriveAxis);

				if (DistanceTo0 > DistanceTo1)
				{
					ensure(!NewDrive.TorqueNodeIndices[1].Contains(i));
					NewDrive.TorqueNodeIndices[1].Add(i);
				}
				else
				{
					ensure(!NewDrive.TorqueNodeIndices[0].Contains(i));
					NewDrive.TorqueNodeIndices[0].Add(i);
				}
			}
		}
	}

	for (const FMHChunkHydraulic& Hydraulic : Chunk.Hydraulics)
	{
		// Look for edge
		const int32 NodeIndices[2] = {
			NodeOffset + Hydraulic.NodeIndices[0],
			NodeOffset + Hydraulic.NodeIndices[1]
		};

		int32 MatchingEdgeIndex = -1;

		for (int32 EdgeIndex = NewMeshInfo.EdgeIndex; EdgeIndex < Edges.Num(); ++EdgeIndex)
		{
			const FMHEdge& Edge = Edges[EdgeIndex];

			if ((Edge.NodeIndices[0] == NodeIndices[0] && Edge.NodeIndices[1] == NodeIndices[1]) ||
				(Edge.NodeIndices[1] == NodeIndices[0] && Edge.NodeIndices[0] == NodeIndices[1]))
			{
				MatchingEdgeIndex = EdgeIndex;
				break;
			}
		}

		if (MatchingEdgeIndex == -1)
		{
			ensure(0);
			continue;
		}

		Hydraulics.AddDefaulted(1);

		FMHHydraulic& NewHydraulic = Hydraulics.Last();
		NewHydraulic.Name = Hydraulic.Name;
		NewHydraulic.NodeIndices[0] = NodeIndices[0];
		NewHydraulic.NodeIndices[1] = NodeIndices[1];
		NewHydraulic.EdgeIndex = MatchingEdgeIndex;
		NewHydraulic.DefaultLength = Edges[MatchingEdgeIndex].DefaultLength;
		NewHydraulic.Length = NewHydraulic.DefaultLength;

		NewMeshInfo.HydraulicNames.Add(Hydraulic.Name);
	}

	NewMeshInfo.NumNodes = Nodes.Num() - NewMeshInfo.NodeIndex;
	NewMeshInfo.NumEdges = Edges.Num() - NewMeshInfo.EdgeIndex;
	NewMeshInfo.NumTriangles = Triangles.Num() - NewMeshInfo.TriangleIndex;
	NewMeshInfo.NumDrives = Drives.Num() - NewMeshInfo.DriveIndex;
	NewMeshInfo.NumHydraulics = Hydraulics.Num() - NewMeshInfo.HydraulicIndex;

	FMHMesh NewMesh;
	NewMesh.Info = NewMeshInfo;

	Meshes.Add(NewMesh);

	return NewMeshInfo;
}

static void _DetectCollisionTriangleToTriangle(const TArray<FMHNode>& Nodes, 
	const FMHTriangle& Triangle1, const FMHTriangle& Triangle2, 
	const int32 TriangleIndex1, const int32 TriangleIndex2,
	TArray<FMHContact>& OutContacts)
{
	SCOPE_CYCLE_COUNTER(STAT_CollisionDetectTriangleToTriangle);

	if (CVarMHPhysicsEnableStaticNodeCollision.GetValueOnAnyThread() == 0)
	{
		if (Nodes[Triangle1.NodeIndices[0]].Mass == 0 &&
			Nodes[Triangle1.NodeIndices[1]].Mass == 0 &&
			Nodes[Triangle1.NodeIndices[2]].Mass == 0)
		{
			return;
		}
	}

	const FVector Positions[3] = {
		Nodes[Triangle1.NodeIndices[0]].Position,
		Nodes[Triangle1.NodeIndices[1]].Position,
		Nodes[Triangle1.NodeIndices[2]].Position
	};

	const FPlane& PrevPlane2 = Triangle2.CachedPrevPlane;
	
	if ((Triangle2.CacheFlags & FMHTriangle::ValidPrevPlane) == 0)
	{
		Triangle2.CachedPrevPlane = FPlane(
			Nodes[Triangle2.NodeIndices[0]].PrevPosition,
			Nodes[Triangle2.NodeIndices[1]].PrevPosition,
			Nodes[Triangle2.NodeIndices[2]].PrevPosition);

		Triangle2.CacheFlags |= FMHTriangle::ValidPrevPlane;
	}

	const FPlane& Plane2 = Triangle2.CachedPlane;
	
	if ((Triangle2.CacheFlags & FMHTriangle::ValidPlane) == 0)
	{
		Triangle2.CachedPlane = FPlane(
			Nodes[Triangle2.NodeIndices[0]].Position,
			Nodes[Triangle2.NodeIndices[1]].Position,
			Nodes[Triangle2.NodeIndices[2]].Position);

		Triangle2.CacheFlags |= FMHTriangle::ValidPlane;
	}

	const float Depth = 5.0f;

	const FVector PrevPositions[3] = {
		Nodes[Triangle1.NodeIndices[0]].PrevPosition + (FVector(PrevPlane2) * Depth),
		Nodes[Triangle1.NodeIndices[1]].PrevPosition + (FVector(PrevPlane2) * Depth),
		Nodes[Triangle1.NodeIndices[2]].PrevPosition + (FVector(PrevPlane2) * Depth)
	};

	// Node to Triangle
	for (int32 i = 0; i < 3; ++i)
	{
		if (Triangle2.HasNodeIndex(Triangle1.NodeIndices[i]))
		{
			// Skip neighboring triangle
			continue;
		}

		const float PrevDistance = PrevPlane2.PlaneDot(PrevPositions[i]);
		const float Distance = Plane2.PlaneDot(Positions[i]);

		if (PrevDistance >= 0 && Distance <= 0)
		{
			float t;

			if (_RayTriangleIntersect(Positions[i], FVector(Plane2),
				Nodes[Triangle2.NodeIndices[0]].Position,
				Nodes[Triangle2.NodeIndices[1]].Position,
				Nodes[Triangle2.NodeIndices[2]].Position, FVector(Plane2), t))
			{
				FMHContact Contact;
				Contact.TriangleIndices[0] = TriangleIndex1;
				Contact.TriangleIndices[1] = TriangleIndex2;
				Contact.Depth = -Distance;
				Contact.Normal = FVector(Plane2);

				Contact.Type = FMHContact::EType::NodeToTriangle;
				Contact.NodeToTriangle.NodeIndex = Triangle1.NodeIndices[i];

				OutContacts.Add(Contact);
			}
		}
	}
}

static void _DetectCollision(
	const TArray<FMHNode>& Nodes, const TArray<FMHTriangle>& Triangles, TArray<FMHMesh>& Meshes,
	TArray<FMHContact>& OutContacts)
{
	SCOPE_CYCLE_COUNTER(STAT_CollisionDetect);

	const int32 NumNodes = Nodes.Num();
	const int32 NumTriangles = Triangles.Num();
	const int32 NumMeshes = Meshes.Num();

	for (int32 NodeIndex = 0; NodeIndex < NumNodes; ++NodeIndex)
	{
		const FMHNode& Node = Nodes[NodeIndex];

		if (Node.Mass == 0.0f)
		{
			continue;
		}

		for (int32 MeshIndex = 0; MeshIndex < NumMeshes; ++MeshIndex)
		{
			const bool bSameMesh = (Node.MeshIndex == MeshIndex);

			if (!Node.bSelfCollision && bSameMesh)
			{
				// Ignore self collision
				continue;
			}

			const FMHMesh& Mesh = Meshes[MeshIndex];
			const int32 NumTriangles = Mesh.Info.NumTriangles;
			const int32 TriangleIndexOffset = Mesh.Info.TriangleIndex;

			for (int32 i = 0; i < NumTriangles; ++i)
			{
				const int32 TriangleIndex = TriangleIndexOffset + i;

				const FMHTriangle& Triangle = Triangles[TriangleIndex];

				if (!Triangle.CachedBBox.Intersect(Node.CachedBBox))
				{
					// Ignore bbox not overlapped
					continue;
				}

				if (bSameMesh && Triangle.HasNodeIndex(NodeIndex))
				{
					// Ignore node within same triangle
					continue;
				}

				SCOPE_CYCLE_COUNTER(STAT_CollisionDetectNodeToTriangle);

				if ((Triangle.CacheFlags & FMHTriangle::ValidPrevPlane) == 0)
				{
					Triangle.CachedPrevPlane = FPlane(
						Nodes[Triangle.NodeIndices[0]].PrevPosition,
						Nodes[Triangle.NodeIndices[1]].PrevPosition,
						Nodes[Triangle.NodeIndices[2]].PrevPosition);

					Triangle.CacheFlags |= FMHTriangle::ValidPrevPlane;
				}

				if ((Triangle.CacheFlags & FMHTriangle::ValidPlane) == 0)
				{
					Triangle.CachedPlane = FPlane(
						Nodes[Triangle.NodeIndices[0]].Position,
						Nodes[Triangle.NodeIndices[1]].Position,
						Nodes[Triangle.NodeIndices[2]].Position);

					Triangle.CacheFlags |= FMHTriangle::ValidPlane;
				}

				const FPlane& PrevPlane = Triangle.CachedPrevPlane;
				const FPlane& Plane = Triangle.CachedPlane;

				const float Depth = 5.0f;
				const FVector PrevPosition = Node.PrevPosition + (FVector(PrevPlane) * Depth);

				const float PrevDistance = PrevPlane.PlaneDot(PrevPosition);
				const float Distance = Plane.PlaneDot(Node.Position);

				if (PrevDistance >= 0 && Distance <= 0)
				{
					float t;

					if (_RayTriangleIntersect(Node.Position, FVector(Plane),
						Nodes[Triangle.NodeIndices[0]].Position,
						Nodes[Triangle.NodeIndices[1]].Position,
						Nodes[Triangle.NodeIndices[2]].Position, FVector(Plane), t))
					{
						FMHContact Contact;
						Contact.TriangleIndices[0] = -1;
						Contact.TriangleIndices[1] = TriangleIndex;
						Contact.Depth = -Distance;
						Contact.Normal = FVector(Plane);

						Contact.Type = FMHContact::EType::NodeToTriangle;
						Contact.NodeToTriangle.NodeIndex = NodeIndex;

						OutContacts.Add(Contact);
					}
				}
			}
		}
	}
}

void FMHPhysics::Tick(float DeltaSeconds)
{
	SCOPE_CYCLE_COUNTER(STAT_MHPhysicsTick);

	Profile.PreTick();

	DeltaSeconds *= CVarMHPhysicsSpeed.GetValueOnAnyThread();

	const float MaxDetalSeconds = 0.03f;
	DeltaSeconds = FMath::Min(MaxDetalSeconds, DeltaSeconds);

	const float fps = FMath::Max(CVarMHPhysicsFPS.GetValueOnAnyThread(), 1.0f);
	const float StepTickSeconds = 1.0f / fps;

	DeltaSeconds += TickSecondLeft;

	while (DeltaSeconds > StepTickSeconds)
	{
		Step(StepTickSeconds);
		DeltaSeconds -= StepTickSeconds;

		++Profile.NumStepsOnLastTick;
	}

	TickSecondLeft = DeltaSeconds;
}

void FMHPhysics::Step(float DeltaSeconds)
{
	SCOPE_CYCLE_COUNTER(STAT_MHPhysicsStep);

	DeltaSeconds = FMath::Clamp(DeltaSeconds, 0.0f, 0.01f);

	Contacts.Reset();

	_DetectCollision(Nodes, Triangles, Meshes, Contacts);

	{
		SCOPE_CYCLE_COUNTER(STAT_Gravity);

		// Apply gravity
		for (FMHNode& Node : Nodes)
		{
			Node.Force = Node.Mass * Setting.Gravity;
		}
	}

	{
		// Apply Hydraulic
		for (const FMHHydraulic& Hydraulic : Hydraulics)
		{
			FMHEdge& Edge = Edges[Hydraulic.EdgeIndex];

			Edge.DefaultLength = Hydraulic.Length;
		}
	}

	{
		SCOPE_CYCLE_COUNTER(STAT_EdgeSpring);

		// Apply spring force
		const int32 NumEdges = Edges.Num();

		struct FHotData
		{
			FVector PositionDiff;
			FVector VelocityDiff;
		};
		TArray<FHotData> HotData;
		HotData.SetNum(NumEdges);

		for (int32 EdgeIndex = 0; EdgeIndex < NumEdges; ++EdgeIndex)
		{
			const FMHEdge& Edge = Edges[EdgeIndex];
			HotData[EdgeIndex].PositionDiff = Nodes[Edge.NodeIndices[0]].Position - Nodes[Edge.NodeIndices[1]].Position;
			HotData[EdgeIndex].VelocityDiff = Nodes[Edge.NodeIndices[0]].Velocity - Nodes[Edge.NodeIndices[1]].Velocity;
		}

		for (int32 EdgeIndex = 0; EdgeIndex < NumEdges; ++EdgeIndex)
		{
			const FMHEdge& Edge = Edges[EdgeIndex];

			const FVector PositionDiff = HotData[EdgeIndex].PositionDiff;
			const FVector VelocityDiff = HotData[EdgeIndex].VelocityDiff;
			const float InvDistance = FMath::InvSqrtEst(PositionDiff.SizeSquared());
			//const FVector Normal = (Distance > SMALL_NUMBER) ? PositionDiff / Distance : FVector::ZeroVector;
			const FVector Normal = PositionDiff * InvDistance;
			const float Distance = PositionDiff | Normal;
			const float Speed = VelocityDiff | Normal;

			const FVector ForceByK = (Distance - Edge.DefaultLength) * Edge.SpringK * Normal;
			const FVector ForceByD = Speed * Edge.SpringD * Normal;

			Nodes[Edge.NodeIndices[0]].Force -= ForceByK + ForceByD;
			Nodes[Edge.NodeIndices[1]].Force += ForceByK + ForceByD;
		}

		// Apply drive force
		for (const FMHDrive& Drive : Drives)
		{
			const FVector DrivePositions[2] = {
				Nodes[Drive.NodeIndices[0]].Position,
				Nodes[Drive.NodeIndices[1]].Position
			};

			const FVector DriveVector = DrivePositions[0] - DrivePositions[1];

			const float DriveTorque = Drive.Torque;
			const float DriveTorque0PerNode = Drive.TorqueNodeIndices[0].Num() > 0 ? 
				DriveTorque / Drive.TorqueNodeIndices[0].Num() : 0.0f;
			const float DriveTorque1PerNode = Drive.TorqueNodeIndices[1].Num() > 0 ? 
				DriveTorque / Drive.TorqueNodeIndices[1].Num() : 0.0f;
			
			for (int32 i : Drive.TorqueNodeIndices[0])
			{
				const FVector& Position = Nodes[i].Position;
				const FVector DriveToNode = Position - DrivePositions[0];
				const FVector ForceVector = (DriveVector ^ DriveToNode).GetSafeNormal();
				const FVector Normal = (ForceVector ^ DriveVector).GetSafeNormal();

				const float Distance = DriveToNode | Normal;
				const float Force = (Distance > SMALL_NUMBER) ? DriveTorque0PerNode / Distance : 0.0f;

				Nodes[i].Force += Force * ForceVector;
			}

			for (int32 i : Drive.TorqueNodeIndices[1])
			{
				const FVector& Position = Nodes[i].Position;
				const FVector DriveToNode = Position - DrivePositions[1];
				const FVector ForceVector = (DriveVector ^ DriveToNode).GetSafeNormal();
				const FVector Normal = (ForceVector ^ DriveVector).GetSafeNormal();

				const float Distance = DriveToNode | Normal;
				const float Force = (Distance > SMALL_NUMBER) ? DriveTorque1PerNode / Distance : 0.0f;

				Nodes[i].Force += Force * ForceVector;
			}
		}
	}

	{
		SCOPE_CYCLE_COUNTER(STAT_CollisionResolve);

		// Add friction forces
		for (const FMHContact& Contact : Contacts)
		{
			if (Contact.Type == FMHContact::EType::NodeToTriangle)
			{
				const float FRICTION_COEFF = CVarMHPhysicsFrictionCoeff.GetValueOnAnyThread(); // TODO: parameterize

				FMHNode& Node = Nodes[Contact.NodeToTriangle.NodeIndex];
				if (Node.Mass != 0.0f)
				{
					const float VertialLoad = FMath::Max(-(Node.Force | Contact.Normal), 0.0f);
					const FVector FrictionDir = ((Node.Velocity ^ Contact.Normal) ^ Contact.Normal).GetSafeNormal();
					const float FrictionStrength = VertialLoad * FRICTION_COEFF;

					Node.Force += FrictionDir * FrictionStrength;
				}

				{
					const int32 TriangleIndex = Contact.TriangleIndices[1];
					const FMHTriangle& Triangle = Triangles[TriangleIndex];
					FMHNode* TriangleNodes[] = {
						&Nodes[Triangle.NodeIndices[0]],
						&Nodes[Triangle.NodeIndices[1]],
						&Nodes[Triangle.NodeIndices[2]]
					};

					for (int32 i = 0; i < 3; ++i)
					{
						FMHNode& TriangleNode = *TriangleNodes[i];

						const float VertialLoad = FMath::Max(TriangleNode.Force | Contact.Normal, 0.0f);
						const FVector FrictionDir = ((TriangleNode.Velocity ^ Contact.Normal) ^ Contact.Normal).GetSafeNormal();
						const float FrictionStrength = VertialLoad * FRICTION_COEFF;

						TriangleNode.Force += FrictionDir * FrictionStrength;
					}
				}
			}
		}

		// Apply Contact impulses
		for (const FMHContact& Contact : Contacts)
		{
			if (Contact.Type == FMHContact::EType::NodeToTriangle)
			{
				FMHNode& Node = Nodes[Contact.NodeToTriangle.NodeIndex];
				if (Node.Mass > SMALL_NUMBER)
				{
					const int32 TriangleIndex = Contact.TriangleIndices[1];
					const FMHTriangle& Triangle = Triangles[TriangleIndex];
					FMHNode* TriangleNodes[] = {
						&Nodes[Triangle.NodeIndices[0]],
						&Nodes[Triangle.NodeIndices[1]],
						&Nodes[Triangle.NodeIndices[2]]
					};

					const float TriangleMass = TriangleNodes[0]->Mass + TriangleNodes[1]->Mass + TriangleNodes[2]->Mass;
					if (TriangleMass > SMALL_NUMBER)
					{
						// Dynamic node to Dynamic triangle
						const float NodeKineticEnergy = (Node.Velocity | Contact.Normal) * Node.Mass;
						const float TriangleKineticEnergy = 
							  (TriangleNodes[0]->Velocity | Contact.Normal) * TriangleNodes[0]->Mass 
							+ (TriangleNodes[1]->Velocity | Contact.Normal) * TriangleNodes[1]->Mass
							+ (TriangleNodes[2]->Velocity | Contact.Normal) * TriangleNodes[2]->Mass;

						const float TotalKineticEnergy = NodeKineticEnergy + TriangleKineticEnergy;
						const float FinalVelocity = TotalKineticEnergy / (Node.Mass + TriangleMass);

						Node.Velocity += (-(Node.Velocity | Contact.Normal) * Contact.Normal) + (FinalVelocity * Contact.Normal);
						TriangleNodes[0]->Velocity += (-(TriangleNodes[0]->Velocity | Contact.Normal) * Contact.Normal) + (FinalVelocity * Contact.Normal);
						TriangleNodes[1]->Velocity += (-(TriangleNodes[1]->Velocity | Contact.Normal) * Contact.Normal) + (FinalVelocity * Contact.Normal);
						TriangleNodes[2]->Velocity += (-(TriangleNodes[2]->Velocity | Contact.Normal) * Contact.Normal) + (FinalVelocity * Contact.Normal);
					}
					else
					{
						// Dynamic node collide to Static Triangle
						Node.Velocity -= FMath::Min(Node.Velocity | Contact.Normal, 0.0f) * Contact.Normal;
					}
				}
			}
		}
		
		// Backward position so that there is no penetration
		for (const FMHContact& Contact : Contacts)
		{
			if (Contact.Type == FMHContact::EType::NodeToTriangle)
			{
				FMHNode& Node = Nodes[Contact.NodeToTriangle.NodeIndex];
				if (Node.Mass > SMALL_NUMBER)
				{
					{
						//Node.Position += Contact.Normal * Contact.Depth;
						const FVector Movement = Node.Position - Node.PrevPosition;
						const float MoveDistance = Movement.Size();
						const float ContactNormalDotMovement = -(Movement.GetSafeNormal() | Contact.Normal);
						const float MoveBackward = ContactNormalDotMovement > SMALL_NUMBER ?
							FMath::Min(MoveDistance, Contact.Depth / ContactNormalDotMovement) : 0.0f;
						
						Node.Position -= Movement.GetSafeNormal() * MoveBackward;
						Node.Force -= FMath::Min(Node.Force | Contact.Normal, 0.0f) * Contact.Normal;
						//Node.Velocity -= FMath::Min(Node.Velocity | Contact.Normal, 0.0f) * Contact.Normal;
					}
				}

				//if (CVarMHPhysicsEnableStaticNodeCollision.GetValueOnAnyThread() != 0)
				{
					const int32 TriangleIndex = Contact.TriangleIndices[1];
					const FMHTriangle& Triangle = Triangles[TriangleIndex];
					FMHNode* TriangleNodes[] = {
						&Nodes[Triangle.NodeIndices[0]],
						&Nodes[Triangle.NodeIndices[1]],
						&Nodes[Triangle.NodeIndices[2]]
					};

					for (int32 i = 0; i < 3; ++i)
					{
						const FVector Movement = TriangleNodes[i]->Position - TriangleNodes[i]->PrevPosition;
						const float MoveDistance = Movement.Size();
						const float ContactNormalDotMovement = Movement.GetSafeNormal() | Contact.Normal;
						const float MoveBackward = ContactNormalDotMovement > SMALL_NUMBER ?
							FMath::Min(MoveDistance, Contact.Depth / ContactNormalDotMovement) : 0.0f;

						TriangleNodes[i]->Position -= Movement.GetSafeNormal() * MoveBackward;
						TriangleNodes[i]->Force -= FMath::Max(TriangleNodes[i]->Force | Contact.Normal, 0.0f) * Contact.Normal;
						//TriangleNodes[i]->Velocity -= FMath::Max(TriangleNodes[i]->Velocity | Contact.Normal, 0.0f) * Contact.Normal;
					}
				}
			}
		}
	}

	// Swap frame
	for (FMHNode& Node : Nodes)
	{
		Node.PrevPosition = Node.Position;
	}

	// Integration
	{
		SCOPE_CYCLE_COUNTER(STAT_Integration);

		const float MaxSpeed = CVarMHPMaxSpeed.GetValueOnAnyThread();
		const float MaxSpeedSquared = FMath::Square(MaxSpeed);

		// Euler Integration
		const int32 NumNodes = Nodes.Num();
		for (int32 NodeIndex = 0; NodeIndex < NumNodes; ++NodeIndex)
		{
			FMHNode& Node = Nodes[NodeIndex];

			if (Node.Mass < SMALL_NUMBER)
			{
				continue;
			}
			FVector Acceleration = Node.Force / Node.Mass;
			FVector AddVelocity = Acceleration * DeltaSeconds;
			Node.Position += (Node.Velocity + (AddVelocity * 0.5f)) * DeltaSeconds;
			Node.Velocity += AddVelocity;

			const float SpeedSquared = Node.Velocity.SizeSquared();
			if (SpeedSquared > MaxSpeedSquared + SMALL_NUMBER)
			{
				Node.Velocity *= MaxSpeed * FMath::InvSqrt(SpeedSquared);

				Profile.TooFastNodes.AddUnique(NodeIndex);
			}
		}
	}

	// Update Caches

	{
		SCOPE_CYCLE_COUNTER(STAT_UpdateNodeBBox);

		for (FMHNode& Node : Nodes)
		{
			Node.UpdateBBox();
		}
	}

	{
		SCOPE_CYCLE_COUNTER(STAT_UpdateTriangleBBox);

		for (FMHTriangle& Triangle : Triangles)
		{
			if (Triangle.CacheFlags & FMHTriangle::ValidPlane)
			{
				Triangle.CachedPrevPlane = Triangle.CachedPlane;
				Triangle.CacheFlags &= ~FMHTriangle::ValidPlane;
				Triangle.CacheFlags |= FMHTriangle::ValidPrevPlane;
			}
			else
			{
				Triangle.CacheFlags &= ~FMHTriangle::ValidPrevPlane;
			}

			Triangle.UpdateBBox(Nodes);
		}
	}
}

const FMHNode* FMHPhysics::FindNode(int32 NodeIndex) const
{
	if (Nodes.IsValidIndex(NodeIndex))
	{
		return &Nodes[NodeIndex];
	}

	return nullptr;
}

const FMHTriangle* FMHPhysics::FindTriangle(int32 TriangleIndex) const
{
	if (Triangles.IsValidIndex(TriangleIndex))
	{
		return &Triangles[TriangleIndex];
	}

	return nullptr;
}

void FMHPhysics::SetDriveTorque(int32 DriveIndex, float Torque)
{
	if (ensure(Drives.IsValidIndex(DriveIndex)))
	{
		Drives[DriveIndex].Torque = Torque;
	}
}

void FMHPhysics::SetHydraulicScale(int32 HydraulicIndex, float Scale)
{
	if (ensure(Hydraulics.IsValidIndex(HydraulicIndex)))
	{
		Hydraulics[HydraulicIndex].Length = Hydraulics[HydraulicIndex].DefaultLength * Scale;
	}
}

void FMHPhysics::DebugDraw(UWorld* World)
{
	SCOPE_CYCLE_COUNTER(STAT_Draw);

	if (CVarMHPhysicsDraw.GetValueOnAnyThread() == 0)
	{
		return;
	}

	GEngine->AddOnScreenDebugMessage((uint64)this, 0, FColor::White,
		FString::Printf(TEXT("[MHPhysics] Nodes: %d, Edges: %d, Triangles: %d"),
			Nodes.Num(), Edges.Num(), Triangles.Num()));

	GEngine->AddOnScreenDebugMessage((uint64)(this + 1), 0, FColor::White,
		FString::Printf(TEXT("[MHPhysics] NumStepsOnLastTick: %03d"),
			Profile.NumStepsOnLastTick));

	for (const FMHNode& Node : Nodes)
	{
		::DrawDebugPoint(World, Node.Position, 5.0f, FColor::White, false, 0);
	}

	for (const FMHEdge& Edge : Edges)
	{
		::DrawDebugLine(World,
			Nodes[Edge.NodeIndices[0]].Position,
			Nodes[Edge.NodeIndices[1]].Position, FColor::White, false, 0);
	}

	for (const FMHTriangle& Triangle : Triangles)
	{
		//::DrawDebugLine(World,
		//	Nodes[Triangle.NodeIndices[0]].Position,
		//	Nodes[Triangle.NodeIndices[1]].Position, FColor::Green, false, 0);

		//::DrawDebugLine(World,
		//	Nodes[Triangle.NodeIndices[0]].Position,
		//	Nodes[Triangle.NodeIndices[1]].Position, FColor::Green, false, 0);

		//::DrawDebugLine(World,
		//	Nodes[Triangle.NodeIndices[0]].Position,
		//	Nodes[Triangle.NodeIndices[2]].Position, FColor::Green, false, 0);

		//::DrawDebugLine(World,
		//	Nodes[Triangle.NodeIndices[1]].Position,
		//	Nodes[Triangle.NodeIndices[2]].Position, FColor::Green, false, 0);

		//FPlane Plane(
		//	Nodes[Triangle.NodeIndices[0]].Position,
		//	Nodes[Triangle.NodeIndices[1]].Position,
		//	Nodes[Triangle.NodeIndices[2]].Position);

		//FVector Center = FVector(
		//	Nodes[Triangle.NodeIndices[0]].Position +
		//	Nodes[Triangle.NodeIndices[1]].Position +
		//	Nodes[Triangle.NodeIndices[2]].Position) * 0.3333333f;

		//::DrawDebugDirectionalArrow(World, Center,
		//	Center + (FVector(Plane) * 30.0f), 5.0f, FColor::White, false, -1.0f, 1, 2.0f);
	}

	for (const FMHDrive& Drive : Drives)
	{
		FVector Location = Nodes[Drive.NodeIndices[0]].Position;
		FString Name = Drive.Name + TEXT("_0");

		::DrawDebugString(World, Nodes[Drive.NodeIndices[0]].Position, Drive.Name + TEXT("_0"), nullptr, FColor::Green, 0.0f);
		::DrawDebugString(World, Nodes[Drive.NodeIndices[1]].Position, Drive.Name + TEXT("_1"), nullptr, FColor::Green, 0.0f);

		::DrawDebugPoint(World, Nodes[Drive.NodeIndices[0]].Position, 10.0f, FColor::Green, false, 0);
		::DrawDebugPoint(World, Nodes[Drive.NodeIndices[1]].Position, 10.0f, FColor::Green, false, 0);

		::DrawDebugLine(World,
			Nodes[Drive.NodeIndices[0]].Position,
			Nodes[Drive.NodeIndices[1]].Position, FColor::Green, false, 0);

		for (int32 i : Drive.TorqueNodeIndices[0])
		{
			::DrawDebugString(World, Nodes[i].Position, Drive.Name + TEXT("_0"), nullptr, FColor::Red, 0.0f);
		}

		for (int32 i : Drive.TorqueNodeIndices[1])
		{
			::DrawDebugString(World, Nodes[i].Position, Drive.Name + TEXT("_1"), nullptr, FColor::Red, 0.0f);
		}
	}

	for (const FMHContact& Contact : Contacts)
	{
		::DrawDebugPoint(World, Nodes[Contact.NodeToTriangle.NodeIndex].Position, 10.0f, FColor::Red, false, 0);
		::DrawDebugDirectionalArrow(World, Nodes[Contact.NodeToTriangle.NodeIndex].Position,
			Nodes[Contact.NodeToTriangle.NodeIndex].Position + (Contact.Normal * 50.0f), 20.0f, FColor::Red, false, -1.0f, 2, 4.0f);
	}

	for (int32 NodeIndex : Profile.TooFastNodes)
	{
		::DrawDebugPoint(World, Nodes[NodeIndex].Position, 15.0f, FColor::Purple, false, 0);
	}
}

#if WITH_EDITOR

// reference : FFbxDataConverter::ConvertPos
FVector _ConvertFbxPos(FbxVector4 Vector)
{
	FVector Out;
	Out[0] = Vector[0];
	// flip Y, then the right-handed axis system is converted to LHS
	Out[1] = -Vector[1];
	Out[2] = Vector[2];
	return Out;
}

static void _FillFbxNodes(FbxNode* Node, FbxNodeAttribute::EType Type, TArray<FbxNode*>& OutNodes)
{
	FbxNodeAttribute::EType FbxNodeAttributeType = Node->GetNodeAttribute() ? 
		Node->GetNodeAttribute()->GetAttributeType() : FbxNodeAttribute::eUnknown;

	if (FbxNodeAttributeType == Type)
	{
		OutNodes.Add(Node);
	}

	for (int32 ChildIndex = 0; ChildIndex < Node->GetChildCount(); ++ChildIndex)
	{
		_FillFbxNodes(Node->GetChild(ChildIndex), Type, OutNodes);
	}
}

// reference: FFbxImporter::ComputeTotalMatrix
static FbxAMatrix _ComputeTotalMatrix(FbxScene* Scene, FbxNode* Node)
{
	FbxAMatrix Geometry;
	FbxVector4 Translation, Rotation, Scaling;
	Translation = Node->GetGeometricTranslation(FbxNode::eSourcePivot);
	Rotation = Node->GetGeometricRotation(FbxNode::eSourcePivot);
	Scaling = Node->GetGeometricScaling(FbxNode::eSourcePivot);
	Geometry.SetT(Translation);
	Geometry.SetR(Rotation);
	Geometry.SetS(Scaling);

	//For Single Matrix situation, obtain transfrom matrix from eDESTINATION_SET, which include pivot offsets and pre/post rotations.
	FbxAMatrix& GlobalTransform = Scene->GetAnimationEvaluator()->GetNodeGlobalTransform(Node);

	const bool bTransformVertexToAbsolute = true;

	//We can bake the pivot only if we don't transform the vertex to the absolute position
	if (!bTransformVertexToAbsolute)
	{
		//if (ImportOptions->bBakePivotInVertex)
		//{
		//	FbxAMatrix PivotGeometry;
		//	FbxVector4 RotationPivot = Node->GetRotationPivot(FbxNode::eSourcePivot);
		//	FbxVector4 FullPivot;
		//	FullPivot[0] = -RotationPivot[0];
		//	FullPivot[1] = -RotationPivot[1];
		//	FullPivot[2] = -RotationPivot[2];
		//	PivotGeometry.SetT(FullPivot);
		//	Geometry = Geometry * PivotGeometry;
		//}
		//else
		{
			//No Vertex transform and no bake pivot, it will be the mesh as-is.
			Geometry.SetIdentity();
		}
	}
	//We must always add the geometric transform. Only Max use the geometric transform which is an offset to the local transform of the node
	FbxAMatrix TotalMatrix = bTransformVertexToAbsolute ? GlobalTransform * Geometry : Geometry;

	return TotalMatrix;
}

bool FMHChunk::LoadFromFbx(const FString& Filename)
{
	Clear();

	UnFbx::FFbxImporter* FFbxImporter = UnFbx::FFbxImporter::GetInstance();

	if (!FFbxImporter->ImportFromFile(*Filename, FPaths::GetExtension(Filename), true))
	{
		return false;
	}

	TArray<FbxNode*> FbxNullNodes;
	_FillFbxNodes(FFbxImporter->Scene->GetRootNode(), FbxNodeAttribute::eNull, FbxNullNodes);

	
	TArray<FbxNode*> FbxMeshArray;
	FFbxImporter->FillFbxMeshArray(FFbxImporter->Scene->GetRootNode(), FbxMeshArray, FFbxImporter);

	for (FbxNode* Node : FbxMeshArray)
	{
		FbxMesh* Mesh = Node->GetMesh();
		check(Mesh);

		// Construct the matrices for the conversion from right handed to left handed system
		FbxAMatrix TotalMatrix = _ComputeTotalMatrix(FFbxImporter->Scene, Node);

		int32 NumVertices = Mesh->GetControlPointsCount();

		for (int32 i = 0; i < NumVertices; ++i)
		{
			FbxVector4 FbxPosition = Mesh->GetControlPoints()[i];
			FbxPosition = TotalMatrix.MultT(FbxPosition);

			Nodes.Add(FMHChunkNode({ _ConvertFbxPos(FbxPosition) }));
		}

		int32 NumEdges = Mesh->GetMeshEdgeCount();

		for (int32 i = 0; i < NumEdges; ++i)
		{
			int32 VertexIndices[2];
			Mesh->GetMeshEdgeVertices(i, VertexIndices[0], VertexIndices[1]);

			Edges.Add(FMHChunkEdge({ VertexIndices[0], VertexIndices[1] }));
		}

		//// Ref : UnFbx::FFbxImporter::BuildStaticMeshFromGeometry
		//if (!Mesh->IsTriangleMesh())
		//{
		//	//UE_LOG(LogFbx, Warning, TEXT("Triangulating static mesh %s"), UTF8_TO_TCHAR(Node->GetName()));

		//	const bool bReplace = true;
		//	FbxNodeAttribute* ConvertedNode = FFbxImporter->GetGeometryConverter()->Triangulate(Mesh, bReplace);

		//	if (ConvertedNode != NULL && ConvertedNode->GetAttributeType() == FbxNodeAttribute::eMesh)
		//	{
		//		Mesh = (fbxsdk::FbxMesh*)ConvertedNode;
		//	}
		//	else
		//	{
		//		//AddTokenizedErrorMessage(FTokenizedMessage::Create(EMessageSeverity::Warning, FText::Format(LOCTEXT("Error_FailedToTriangulate", "Unable to triangulate mesh '{0}'"), FText::FromString(Mesh->GetName()))), FFbxErrors::Generic_Mesh_TriangulationFailed);
		//		ensure(0);
		//		return false; // not clean, missing some dealloc
		//	}
		//}

		//if (ensure(Mesh->IsTriangleMesh()))
		{
			int32 NumTriangles = Mesh->GetPolygonCount();

			for (int32 i = 0; i < NumTriangles; ++i)
			{
				if (Mesh->GetPolygonSize(i) == 3)
				{
					bool bValidTriangle = true;
					int32 NodeIndices[3];

					for (int32 CornerIndex = 0; CornerIndex < 3; ++CornerIndex)
					{
						int32 ControlPointIndex = Mesh->GetPolygonVertex(i, CornerIndex);

						if (!ensure(Nodes.IsValidIndex(ControlPointIndex)))
						{
							bValidTriangle = false;
							break;
						}

						NodeIndices[CornerIndex] = ControlPointIndex;
					}

					if (bValidTriangle)
					{
						Triangles.Add(FMHChunkTriangle({ NodeIndices[0], NodeIndices[2], NodeIndices[1] }));
					}
				}
			}
		}
	}

	const static float NODE_SEARCH_DISTANCE_SQUARED = FMath::Square(0.01f);

	struct FNullNode
	{
		FString Name;
		int32 NodeIndex;
	};
	TArray<FNullNode> NullNodes;

	for (FbxNode* NullNode : FbxNullNodes)
	{
		FbxAMatrix TotalMatrix = _ComputeTotalMatrix(FFbxImporter->Scene, NullNode);
		FVector Location = _ConvertFbxPos(TotalMatrix.MultT(FbxVector4(0, 0, 0)));

		for (int32 NodeIndex = 0; NodeIndex < Nodes.Num(); ++NodeIndex)
		{
			const FMHChunkNode& ChunkNode = Nodes[NodeIndex];

			if ((ChunkNode.Position - Location).SizeSquared() < NODE_SEARCH_DISTANCE_SQUARED)
			{
				FString Name = NullNode->GetName();
				NullNodes.Add(FNullNode({ Name, NodeIndex }));
				break;
			}
		}
	}

	// Look for Drives
	for (const FNullNode& NullNode : NullNodes)
	{
		if (NullNode.Name.StartsWith(TEXT("__DRIVE0_")))
		{
			FString PairName = NullNode.Name.Replace(TEXT("__DRIVE0_"), TEXT("__DRIVE1_"));
			FNullNode* PairNode = NullNodes.FindByPredicate([PairName](const FNullNode& Node) {
				return Node.Name == PairName;
			});

			if (ensure(PairNode))
			{
				FMHChunkDrive Drive;
				Drive.Name = NullNode.Name.Replace(TEXT("__DRIVE0_"), TEXT(""));
				Drive.NodeIndices[0] = NullNode.NodeIndex;
				Drive.NodeIndices[1] = PairNode->NodeIndex;
				Drives.Add(Drive);
			}
		}
	}

	// Look for Hydraulics
	for (const FNullNode& NullNode : NullNodes)
	{
		if (NullNode.Name.StartsWith(TEXT("__HYD0_")))
		{
			FString PairName = NullNode.Name.Replace(TEXT("__HYD0_"), TEXT("__HYD1_"));
			FNullNode* PairNode = NullNodes.FindByPredicate([PairName](const FNullNode& Node) {
				return Node.Name == PairName;
			});

			if (ensure(PairNode))
			{
				FMHChunkHydraulic Hydraulic;
				Hydraulic.Name = NullNode.Name.Replace(TEXT("__HYD0_"), TEXT(""));
				Hydraulic.NodeIndices[0] = NullNode.NodeIndex;
				Hydraulic.NodeIndices[1] = PairNode->NodeIndex;
				Hydraulics.Add(Hydraulic);
			}
		}
	}

	return true;
}

#endif // WITH_EDITORONLY_DATA

void FMHTriangle::UpdateBBox(const TArray<FMHNode>& Nodes)
{
	CachedBBox = Nodes[NodeIndices[0]].CachedBBox;
	CachedBBox += Nodes[NodeIndices[1]].CachedBBox;
	CachedBBox += Nodes[NodeIndices[2]].CachedBBox;
}

