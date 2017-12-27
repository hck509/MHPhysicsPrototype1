// Fill out your copyright notice in the Description page of Project Settings.

#include "MHStaticMeshComponent.h"
#include "CustomMeshComponent.h"
#include "MHPhysicsPrototype1GameModeBase.h"
#include "Engine/World.h"

#if WITH_EDITORONLY_DATA

#include "FeedbackContext.h"
#include "DetailLayoutBuilder.h"
#include "DetailCategoryBuilder.h"
#include "IDetailsView.h"
#include "DetailWidgetRow.h"

#include "IDesktopPlatform.h"
#include "DesktopPlatformModule.h"
#include "SlateApplication.h"
#include "EditorDirectories.h"

#define LOCTEXT_NAMESPACE "MHStaticMeshComponentDetails"

TSharedRef<IDetailCustomization> FMHStaticMeshComponentDetails::MakeInstance()
{
	return MakeShareable(new FMHStaticMeshComponentDetails);
}

void FMHStaticMeshComponentDetails::CustomizeDetails(IDetailLayoutBuilder& DetailLayout)
{
	const TArray< TWeakObjectPtr<UObject> >& SelectedObjects = DetailLayout.GetDetailsView()->GetSelectedObjects();

	for (int32 ObjectIndex = 0; ObjectIndex < SelectedObjects.Num(); ++ObjectIndex)
	{
		const TWeakObjectPtr<UObject>& CurrentObject = SelectedObjects[ObjectIndex];
		if (CurrentObject.IsValid())
		{
			UMHStaticMeshComponent* CurrentComponents = Cast<UMHStaticMeshComponent>(CurrentObject.Get());
			if (CurrentComponents != NULL)
			{
				Components.Add(CurrentComponents);
			}

			AActor* Actor = Cast<AActor>(CurrentObject.Get());
			if (Actor)
			{
				TInlineComponentArray<UMHStaticMeshComponent*> Components;
				Actor->GetComponents(Components);

				for (UMHStaticMeshComponent* Component : Components)
				{
					Components.Add(Component);
				}
			}
		}
	}

	DetailLayout.EditCategory("MHPhysics")
		.AddCustomRow(NSLOCTEXT("MHStaticMeshComponentDetails", "Import FBX", "Import FBX"))
		.NameContent()
		[
			SNullWidget::NullWidget
		]
	.ValueContent()
		[
			SNew(SBox)
			.WidthOverride(125)
		[
			SNew(SButton)
			.ContentPadding(3)
		.VAlign(VAlign_Center)
		.HAlign(HAlign_Center)
		.OnClicked(this, &FMHStaticMeshComponentDetails::ImportFBX)
		[
			SNew(STextBlock)
			.Text(NSLOCTEXT("CombatCameraComponentDetails", "Import FBX", "Import FBX"))
		.Font(IDetailLayoutBuilder::GetDetailFont())
		]
		]
		];
}

FReply FMHStaticMeshComponentDetails::ImportFBX()
{
	for (auto& Component : Components)
	{
		if (Component.IsValid())
		{
			Component->ImportFBX();
		}
	}

	return FReply::Handled();
}

#endif // WITH_EDITORONLY_DATA


// Sets default values for this component's properties
UMHStaticMeshComponent::UMHStaticMeshComponent()
{
	// Set this component to be initialized when the game starts, and to be ticked every frame.  You can turn these features
	// off to improve performance if you don't need them.
	PrimaryComponentTick.bCanEverTick = true;

	CustomMeshComponent = CreateDefaultSubobject<UCustomMeshComponent>(TEXT("CustomMeshComponent"));
	CustomMeshComponent->SetupAttachment(this);

	MassInKg = 100.0f;
	SpringK = 1000.0f;
	SpringD = 10.0f;
}

// Called when the game starts
void UMHStaticMeshComponent::BeginPlay()
{
	Super::BeginPlay();
}

void UMHStaticMeshComponent::OnRegister()
{
	Super::OnRegister();

	InitializeFromStaticMesh();
	InitializeFromChunk();
}

// Called every frame
void UMHStaticMeshComponent::TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction)
{
	Super::TickComponent(DeltaTime, TickType, ThisTickFunction);

	UpdateCustomMeshFromMHPhysics();
}

void UMHStaticMeshComponent::InitializeFromStaticMesh()
{
	if (StaticMesh)
	{
		if (ensure(StaticMesh->RenderData.IsValid()) && ensure(StaticMesh->RenderData->LODResources.Num() > 0))
		{
			TArray<FCustomMeshTriangle> Triangles;

			const FStaticMeshLODResources& Resource = StaticMesh->RenderData->LODResources[0];

			FIndexArrayView IndexArrayView = Resource.IndexBuffer.GetArrayView();

			ensure(IndexArrayView.Num() % 3 == 0);

			for (int32 Index = 0; Index + 2 < IndexArrayView.Num(); Index += 3)
			{
				Triangles.Add(FCustomMeshTriangle({
					Resource.PositionVertexBuffer.VertexPosition(IndexArrayView[Index]),
					Resource.PositionVertexBuffer.VertexPosition(IndexArrayView[Index + 1]),
					Resource.PositionVertexBuffer.VertexPosition(IndexArrayView[Index + 2]) }));
			}

			CustomMeshComponent->SetCustomMeshTriangles(Triangles);

			AMHPhysicsPrototype1GameModeBase* GameMode = Cast<AMHPhysicsPrototype1GameModeBase>(GetWorld()->GetAuthGameMode());

			if (GameMode)
			{
				MHMeshInfo = GameMode->GetMHPhyscis().GenerateFromStaticMesh(*StaticMesh, GetComponentTransform(), MassInKg, SpringK, SpringD);
			}
		}
	}
}

void UMHStaticMeshComponent::InitializeFromChunk()
{
	if (MHChunk.Nodes.Num() == 0)
	{
		return;
	}

	UpdateCustomMeshFromChunk();

	AMHPhysicsPrototype1GameModeBase* GameMode = Cast<AMHPhysicsPrototype1GameModeBase>(GetWorld()->GetAuthGameMode());

	if (GameMode)
	{
		MHMeshInfo = GameMode->GetMHPhyscis().GenerateFromChunk(MHChunk, GetComponentTransform(), MassInKg, SpringK, SpringD);
	}
}

void UMHStaticMeshComponent::UpdateCustomMeshFromChunk()
{
	if (MHChunk.Nodes.Num() == 0)
	{
		return;
	}

	TArray<FCustomMeshTriangle> Triangles;

	const FTransform ComponentTransform = GetComponentTransform();

	for (int32 i = 0; i < MHChunk.Triangles.Num(); ++i)
	{
		Triangles.Add(FCustomMeshTriangle({
			MHChunk.Nodes[MHChunk.Triangles[i].NodeIndices[0]].Position,
			MHChunk.Nodes[MHChunk.Triangles[i].NodeIndices[2]].Position,
			MHChunk.Nodes[MHChunk.Triangles[i].NodeIndices[1]].Position }));
	}

	CustomMeshComponent->SetCustomMeshTriangles(Triangles);
}

void UMHStaticMeshComponent::UpdateCustomMeshFromMHPhysics()
{
	if (MHMeshInfo.NumTriangles == 0)
	{
		return;
	}

	AMHPhysicsPrototype1GameModeBase* GameMode = Cast<AMHPhysicsPrototype1GameModeBase>(GetWorld()->GetAuthGameMode());

	if (!GameMode)
	{
		return;
	}

	const FTransform ComponentTransform = GetComponentTransform();

	FMHPhysics& MHPhysics = GameMode->GetMHPhyscis();

	TArray<FCustomMeshTriangle> Triangles;

	for (int32 i = 0; i < MHMeshInfo.NumTriangles; ++i)
	{
		const FMHTriangle* Triangle = MHPhysics.FindTriangle(MHMeshInfo.TriangleIndex + i);

		if (ensure(Triangle))
		{
			const FMHNode* Node0 = MHPhysics.FindNode(Triangle->NodeIndices[0]);
			const FMHNode* Node1 = MHPhysics.FindNode(Triangle->NodeIndices[1]);
			const FMHNode* Node2 = MHPhysics.FindNode(Triangle->NodeIndices[2]);

			if (ensure(Node0 && Node1 && Node2))
			{
				Triangles.Add(FCustomMeshTriangle({
					ComponentTransform.InverseTransformPosition(Node0->Position),
					ComponentTransform.InverseTransformPosition(Node2->Position),
					ComponentTransform.InverseTransformPosition(Node1->Position) }));
			}
		}
	}

	CustomMeshComponent->SetCustomMeshTriangles(Triangles);
}

#if WITH_EDITORONLY_DATA

void UMHStaticMeshComponent::ImportFBX()
{
	TArray<FString> OpenFilenames;
	IDesktopPlatform* DesktopPlatform = FDesktopPlatformModule::Get();
	bool bOpen = false;
	if (DesktopPlatform)
	{
		FString ExtensionStr;
		ExtensionStr += TEXT("FBX (*.fbx)|*.fbx|");

		bOpen = DesktopPlatform->OpenFileDialog(
			FSlateApplication::Get().FindBestParentWindowHandleForDialogs(nullptr),
			NSLOCTEXT("MHStaticMeshComponent", "ImportFBX", "Import FBX from...").ToString(),
			FEditorDirectories::Get().GetLastDirectory(ELastDirectory::FBX),
			TEXT(""),
			*ExtensionStr,
			EFileDialogFlags::None,
			OpenFilenames
		);

		if (bOpen && OpenFilenames.Num() > 0)
		{
			MHChunk.LoadFromFbx(OpenFilenames[0]);
		}
	}
}

#endif // WITH_EDITORONLY_DATA