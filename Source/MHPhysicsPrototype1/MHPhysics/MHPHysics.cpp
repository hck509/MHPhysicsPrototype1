#include "MHPhysics.h"
#include "DrawDebugHelpers.h"
#include "EngineUtils.h"
#include "Engine/StaticMeshActor.h"

static TAutoConsoleVariable<float> CVarMHPhysicsSpeed(
	TEXT("mhp.speed"),
	1.0f,
	TEXT("Speed multiplier")
);

static bool _RayTriangleIntersect(
	const FVector &orig, const FVector &dir,
	const FVector &v0, const FVector &v1, const FVector &v2,
	float &t)
{
	// compute plane's normal
	FVector v0v1 = v1 - v0;
	FVector v0v2 = v2 - v0;
	// no need to normalize
	FVector N = v0v1 ^ v0v2; // N 
	N.Normalize();
	//float area2 = N.Size();

	// Step 1: finding P

	// check if ray and plane are parallel ?
	float NdotRayDirection = N | dir;
	if (FMath::Abs(NdotRayDirection) < SMALL_NUMBER) // almost 0 
		return false; // they are parallel so they don't intersect ! 

	// compute d parameter using equation 2
	float d = N | v0;

	// compute t (equation 3)
	t = (d - (N | orig)) / NdotRayDirection;
	// check if the triangle is in behind the ray
	if (t < 0) return false; // the triangle is behind 

	// compute the intersection point using equation 1
	FVector P = orig + t * dir;

	// Step 2: inside-outside test
	FVector C; // vector perpendicular to triangle's plane 

	// edge 0
	FVector edge0 = v1 - v0;
	FVector vp0 = P - v0;
	C = edge0 ^ vp0;
	if ((N | C) < 0) return false; // P is on the right side 

	// edge 1
	FVector edge1 = v2 - v1;
	FVector vp1 = P - v1;
	C = edge1 ^ vp1;
	if ((N | C) < 0)  return false; // P is on the right side 

	// edge 2
	FVector edge2 = v0 - v2;
	FVector vp2 = P - v2;
	C = edge2 ^ vp2;
	if ((N | C) < 0) return false; // P is on the right side; 

	return true; // this ray hits the triangle 
}


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
		FMHTriangle NewTriangle = FMHTriangle({
			NodeOffset + static_cast<int32>(IndexArrayView[Index]),
			NodeOffset + static_cast<int32>(IndexArrayView[Index + 2]),
			NodeOffset + static_cast<int32>(IndexArrayView[Index + 1]) });

		Triangles.Add(NewTriangle);

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

	const FPlane Plane1(
		Nodes[Triangle1.NodeIndices[0]].Position,
		Nodes[Triangle1.NodeIndices[1]].Position,
		Nodes[Triangle1.NodeIndices[2]].Position);

	const FPlane PrevPlane2(
		Nodes[Triangle2.NodeIndices[0]].PrevPosition,
		Nodes[Triangle2.NodeIndices[1]].PrevPosition,
		Nodes[Triangle2.NodeIndices[2]].PrevPosition);

	const FPlane Plane2(
		Nodes[Triangle2.NodeIndices[0]].Position,
		Nodes[Triangle2.NodeIndices[1]].Position,
		Nodes[Triangle2.NodeIndices[2]].Position);

	const float Depth = 5.0f;

	const FVector PrevPositions[3] = {
		Nodes[Triangle1.NodeIndices[0]].PrevPosition + (FVector(PrevPlane2) * Depth),
		Nodes[Triangle1.NodeIndices[1]].PrevPosition + (FVector(PrevPlane2) * Depth),
		Nodes[Triangle1.NodeIndices[2]].PrevPosition + (FVector(PrevPlane2) * Depth)
	};

	// Node to Plane
	for (int32 i = 0; i < 3; ++i)
	{
		if (Triangle2.HasNodeIndex(Triangle1.NodeIndices[i]))
		{
			continue;
		}

		const float PrevDistance = PrevPlane2.PlaneDot(PrevPositions[i]);
		const float Distance = Plane2.PlaneDot(Positions[i]);

		if (PrevDistance >= 0 && Distance <= 0)
		{
			float t;

			//if (_RayTriangleIntersect(PrevPositions[i], (Positions[i] - PrevPositions[i]).GetSafeNormal(),

			if (_RayTriangleIntersect(Positions[i], FVector(Plane2.X, Plane2.Y, Plane2.Z),
				Nodes[Triangle2.NodeIndices[0]].Position,
				Nodes[Triangle2.NodeIndices[1]].Position,
				Nodes[Triangle2.NodeIndices[2]].Position, t))
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
	DeltaSeconds = FMath::Clamp(DeltaSeconds, 0.0f, 0.01f);
	DeltaSeconds *= CVarMHPhysicsSpeed.GetValueOnAnyThread();

	Contacts.Empty();

	_DetectCollision(Nodes, Triangles, Contacts);

	// Apply gravity
	for (FMHNode& Node : Nodes)
	{
		Node.Force = Node.Mass * Setting.Gravity;
	}

	// Apply spring force
	for (FMHEdge& Edge : Edges)
	{
		const FVector Positions[2] = {
			Nodes[Edge.NodeIndices[0]].Position,
			Nodes[Edge.NodeIndices[1]].Position
		};

		const FVector PositionDiff = Nodes[Edge.NodeIndices[0]].Position - Nodes[Edge.NodeIndices[1]].Position;
		const FVector VelocityDiff = Nodes[Edge.NodeIndices[0]].Velocity - Nodes[Edge.NodeIndices[1]].Velocity;
		const float Distance = PositionDiff.Size();
		const FVector Normal = (Distance > SMALL_NUMBER) ? PositionDiff / Distance : FVector::ZeroVector;
		const float Speed = VelocityDiff | Normal;

		Nodes[Edge.NodeIndices[0]].Force -= (Distance - Edge.DefaultLength) * Edge.SpringK * Normal;
		Nodes[Edge.NodeIndices[0]].Force -= Speed * Edge.SpringD * Normal;
		Nodes[Edge.NodeIndices[1]].Force += (Distance - Edge.DefaultLength) * Edge.SpringK * Normal;
		Nodes[Edge.NodeIndices[1]].Force += Speed * Edge.SpringD * Normal;
	}

	// Adjust position by contact
	for (const FMHContact& Contact : Contacts)
	{
		if (Contact.Type == FMHContact::EType::NodeToTriangle)
		{
			FMHNode& Node = Nodes[Contact.NodeToTriangle.NodeIndex];
			if (Node.Mass != 0.0f)
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
					Node.Velocity -= FMath::Min(Node.Velocity | Contact.Normal, 0.0f) * Contact.Normal;
				}
			}

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
				TriangleNodes[i]->Position -= Movement;
				TriangleNodes[i]->Force -= FMath::Max(TriangleNodes[i]->Force | Contact.Normal, 0.0f) * Contact.Normal;
				TriangleNodes[i]->Velocity -= FMath::Max(TriangleNodes[i]->Velocity | Contact.Normal, 0.0f) * Contact.Normal;
			}
		}
	}

	// Swap frame
	for (FMHNode& Node : Nodes)
	{
		Node.PrevPosition = Node.Position;
	}

	// Euler Integration
	for (FMHNode& Node : Nodes)
	{
		FVector Acceleration = Node.Mass > SMALL_NUMBER ? Node.Force / Node.Mass : FVector::ZeroVector;
		FVector AddVelocity = Acceleration * DeltaSeconds;
		Node.Position += (Node.Velocity + (AddVelocity * 0.5f)) * DeltaSeconds;
		Node.Velocity += AddVelocity;
	}
}

void FMHPhysics::DebugDraw(UWorld* World)
{
	for (const FMHNode& Node : Nodes)
	{
		::DrawDebugPoint(World, Node.Position, 10.0f, FColor::White, false, 0);
	}

	for (const FMHEdge& Edge : Edges)
	{
		::DrawDebugLine(World,
			Nodes[Edge.NodeIndices[0]].Position,
			Nodes[Edge.NodeIndices[1]].Position, FColor::White, false, 0);
	}

	for (const FMHContact& Contact : Contacts)
	{
		::DrawDebugPoint(World, Nodes[Contact.NodeToTriangle.NodeIndex].Position, 10.0f, FColor::Red, false, 0);
		::DrawDebugDirectionalArrow(World, Nodes[Contact.NodeToTriangle.NodeIndex].Position,
			Nodes[Contact.NodeToTriangle.NodeIndex].Position + (Contact.Normal * 30.0f), 5.0f, FColor::Red, false, -1.0f, 1, 2.0f);
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

		FPlane Plane(
			Nodes[Triangle.NodeIndices[0]].Position,
			Nodes[Triangle.NodeIndices[1]].Position,
			Nodes[Triangle.NodeIndices[2]].Position);

		FVector Center = FVector(
			Nodes[Triangle.NodeIndices[0]].Position +
			Nodes[Triangle.NodeIndices[1]].Position +
			Nodes[Triangle.NodeIndices[2]].Position) * 0.3333333f;

		::DrawDebugDirectionalArrow(World, Center,
			Center + (FVector(Plane) * 30.0f), 5.0f, FColor::White, false, -1.0f, 1, 2.0f);
	}
}
