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
			GenerateFromStaticMesh(*StaticMesh, StaticMeshActor->GetActorTransform(), 0.0f);
		}
	}
}

FMHMeshInfo FMHPhysics::GenerateFromStaticMesh(const UStaticMesh& Mesh, const FTransform& Transform, float MeshMassInKg)
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
	const uint32 NumVertices = Resource.PositionVertexBuffer.GetNumVertices();
	const float MassPerNode = NumVertices > 0 ? MeshMassInKg / NumVertices : 0;
	
	for (uint32 VertexIndex = 0; VertexIndex < NumVertices; ++VertexIndex)
	{
		const FVector Position = Transform.TransformPosition(Resource.PositionVertexBuffer.VertexPosition(VertexIndex));

		FMHNode NewNode;
		NewNode.InitNode(MassPerNode, Position);

		Nodes.Add(NewNode);
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

static void _DetectCollisionTriangleToTriangle(const TArray<FMHNode>& Nodes, 
	const FMHTriangle& Triangle1, const FMHTriangle& Triangle2, 
	const int32 TriangleIndex1, const int32 TriangleIndex2,
	TArray<FMHContact>& OutContacts)
{
	const FVector Positions[3] = {
		Nodes[Triangle1.NodeIndices[0]].Position,
		Nodes[Triangle1.NodeIndices[1]].Position,
		Nodes[Triangle1.NodeIndices[2]].Position
	};

	const FVector PrevPositions[3] = {
		Nodes[Triangle1.NodeIndices[0]].PrevPosition,
		Nodes[Triangle1.NodeIndices[1]].PrevPosition,
		Nodes[Triangle1.NodeIndices[2]].PrevPosition
	};

	const FPlane Plane2(
		Nodes[Triangle2.NodeIndices[0]].Position,
		Nodes[Triangle2.NodeIndices[1]].Position,
		Nodes[Triangle2.NodeIndices[2]].Position);

	for (int32 i = 0; i < 3; ++i)
	{
		const float PrevDistance = Plane2.PlaneDot(PrevPositions[i]);
		const float Distance = Plane2.PlaneDot(Positions[i]);

		if (PrevDistance >= 0 && Distance <= 0)
		{
			FMHContact Contact;
			Contact.TriangleIndices[0] = TriangleIndex1;
			Contact.TriangleIndices[1] = TriangleIndex2;
			Contact.Depth = -Distance;
			Contact.Normal = FVector(Plane2.X, Plane2.Y, Plane2.Z);

			Contact.NodeToPlane.NodeIndex = Triangle1.NodeIndices[i];

			OutContacts.Add(Contact);
		}
	}
}

static void _DetectCollision(const TArray<FMHNode>& Nodes, const TArray<FMHTriangle>& Triangles, 
	TArray<FMHContact>& OutContacts)
{
	const int32 NumTriangles = Triangles.Num();

	for (int TriangleIndex1 = 0; TriangleIndex1 < NumTriangles; ++TriangleIndex1)
	{
		for (int TriangleIndex2 = TriangleIndex1 + 1; TriangleIndex2 < NumTriangles; ++TriangleIndex2)
		{
			_DetectCollisionTriangleToTriangle(Nodes, Triangles[TriangleIndex1], Triangles[TriangleIndex2], TriangleIndex1, TriangleIndex2, OutContacts);
			_DetectCollisionTriangleToTriangle(Nodes, Triangles[TriangleIndex2], Triangles[TriangleIndex1], TriangleIndex2, TriangleIndex1, OutContacts);
		}
	}
}

void FMHPhysics::Tick(float DeltaSeconds)
{
	Contacts.Empty();

	_DetectCollision(Nodes, Triangles, Contacts);

	// Adjust position by contact
	for (const FMHContact& Contact : Contacts)
	{
		FMHNode& Node = Nodes[Contact.NodeToPlane.NodeIndex];
		if (Node.Mass != 0.0f)
		{
			Node.Position += Contact.Normal * Contact.Depth;
			Node.Force -= (Node.Force | Contact.Normal) * Contact.Normal;
			Node.Velocity -= (Node.Velocity | Contact.Normal) * Contact.Normal;
		}
	}

	// Apply gravity
	for (FMHNode& Node : Nodes)
	{
		Node.Force = Node.Mass * Setting.Gravity;
	}

	// Swap frame
	for (FMHNode& Node : Nodes)
	{
		Node.PrevPosition = Node.Position;
	}

	// Euler Integration
	for (FMHNode& Node : Nodes)
	{
		FVector Acceleration = Node.Mass > 0.0f ? Node.Force / Node.Mass : FVector::ZeroVector;
		FVector AddVelocity = Acceleration * DeltaSeconds;
		Node.Position += (Node.Velocity + (AddVelocity * 0.5f)) * DeltaSeconds;
		Node.Velocity += AddVelocity;
	}
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
