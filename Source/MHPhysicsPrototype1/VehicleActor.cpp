#include "VehicleActor.h"
#include "Engine.h"
#include "MHPhysicsPrototype1GameModeBase.h"
#include "MHPhysics/MHPHysics.h"

// Sets default values
AVehicleActor::AVehicleActor()
{
 	// Set this pawn to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;

	SceneComponent = CreateDefaultSubobject<USceneComponent>(TEXT("SceneComponent"));
	RootComponent = SceneComponent;

	MeshComponent = CreateDefaultSubobject<UMHStaticMeshComponent>(TEXT("MeshComponent"));
	MeshComponent->SetupAttachment(RootComponent);

	// Create a spring arm component
	CameraSpringArm = CreateDefaultSubobject<USpringArmComponent>(TEXT("CameraSpringArm"));
	CameraSpringArm->TargetOffset = FVector(0.f, 0.f, 100.f);
	CameraSpringArm->SetRelativeRotation(FRotator(-15.f, 0.f, 0.f));
	CameraSpringArm->SetupAttachment(RootComponent);
	CameraSpringArm->TargetArmLength = 800.0f;
	CameraSpringArm->bEnableCameraRotationLag = true;
	CameraSpringArm->CameraRotationLagSpeed = 4.f;
	CameraSpringArm->bInheritPitch = false;
	CameraSpringArm->bInheritRoll = false;

	// Create camera component 
	Camera = CreateDefaultSubobject<UCameraComponent>(TEXT("Camera"));
	Camera->SetupAttachment(CameraSpringArm, USpringArmComponent::SocketName);
	Camera->bUsePawnControlRotation = false;
	Camera->FieldOfView = 90.f;

	Control = CreateDefaultSubobject<UVehicleControl>(TEXT("Control"));

	CenterNodeIndex = -1;
}

// Called when the game starts or when spawned
void AVehicleActor::BeginPlay()
{
	Super::BeginPlay();

	AMHPhysicsPrototype1GameModeBase* GameMode = Cast<AMHPhysicsPrototype1GameModeBase>(GetWorld()->GetAuthGameMode());

	if (GameMode)
	{
		FMHPhysics& Physics = GameMode->GetMHPhyscis();

		// find nearest to center node index
		const FVector CenterPosition = GetActorLocation();
		float NearestDistanceSquared = 10E10;
		int32 NearestIndex = -1;

		const FMHMeshInfo& MeshInfo = MeshComponent->GetMHMeshInfo();
		for (int32 NodeIndex = 0; NodeIndex < MeshInfo.NumNodes; ++NodeIndex)
		{
			const FMHNode* Node = Physics.FindNode(MeshInfo.NodeIndex + NodeIndex);
			if (Node)
			{
				const float DistanceSquared = (Node->Position - CenterPosition).SizeSquared();
				if (NearestDistanceSquared > DistanceSquared)
				{
					NearestDistanceSquared = DistanceSquared;
					NearestIndex = NodeIndex;
				}
			}
		}

		CenterNodeIndex = NearestIndex;
	}
}

// Called every frame
void AVehicleActor::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

	UpdateLocationFromPhysics();

	Control->Tick(DeltaTime);

	UpdateDriveTorque();
	UpdateSteeringHydraulics();

	TickCameraRotation();
}

void AVehicleActor::UpdateLocationFromPhysics()
{
	if (CenterNodeIndex != -1)
	{
		AMHPhysicsPrototype1GameModeBase* GameMode = Cast<AMHPhysicsPrototype1GameModeBase>(GetWorld()->GetAuthGameMode());

		if (GameMode)
		{
			FMHPhysics& Physics = GameMode->GetMHPhyscis();
			const FMHMeshInfo& MeshInfo = MeshComponent->GetMHMeshInfo();

			if (CenterNodeIndex < MeshInfo.NumNodes)
			{
				const FMHNode* Node = Physics.FindNode(MeshInfo.NodeIndex + CenterNodeIndex);
				if (ensure(Node))
				{
					SetActorLocation(Node->Position);
				}
			}
		}
	}
}

void AVehicleActor::UpdateDriveTorque()
{
	AMHPhysicsPrototype1GameModeBase* GameMode = Cast<AMHPhysicsPrototype1GameModeBase>(GetWorld()->GetAuthGameMode());

	if (ensure(GameMode))
	{
		FMHPhysics& MHPhysics = GameMode->GetMHPhyscis();

		for (int32 i = 0; i < MeshComponent->GetMHMeshInfo().NumDrives; ++i)
		{
			const int32 DriveIndex = MeshComponent->GetMHMeshInfo().DriveIndex + i;

			MHPhysics.SetDriveTorque(DriveIndex, 10000000 * Control->GetThrottle());
		}
	}
}

void AVehicleActor::UpdateSteeringHydraulics()
{
	AMHPhysicsPrototype1GameModeBase* GameMode = Cast<AMHPhysicsPrototype1GameModeBase>(GetWorld()->GetAuthGameMode());

	if (ensure(GameMode))
	{
		FMHPhysics& MHPhysics = GameMode->GetMHPhyscis();

		for (int32 i = 0; i < MeshComponent->GetMHMeshInfo().NumHydraulics; ++i)
		{
			const int32 HydraulicIndex = MeshComponent->GetMHMeshInfo().HydraulicIndex + i;

			MHPhysics.SetHydraulicScale(HydraulicIndex, 1.0f + (Control->GetSteer() * 0.2f));
		}
	}
}

void AVehicleActor::TickCameraRotation()
{
	APlayerController* PossessedPlayerController = Cast<APlayerController>(Controller);

	if (PossessedPlayerController && PossessedPlayerController->IsInputKeyDown(EKeys::RightMouseButton))
	{
		float DeltaX, DeltaY;
		PossessedPlayerController->GetInputMouseDelta(DeltaX, DeltaY);

		//CameraSpringArm->AddRelativeRotation(FRotator(0, DeltaX, 0));

		FRotator RotationDelta(0, DeltaX, 0);
		FTransform Transform = CameraSpringArm->GetComponentTransform();
		Transform.SetRotation(RotationDelta.Quaternion() * Transform.GetRotation());
		Transform.NormalizeRotation();
		CameraSpringArm->SetWorldTransform(Transform);
	}
}

void AVehicleActor::PossessedBy(AController* NewController)
{
	Super::PossessedBy(NewController);

	Camera->Activate();
}

// Called to bind functionality to input
void AVehicleActor::SetupPlayerInputComponent(UInputComponent* PlayerInputComponent)
{
	Super::SetupPlayerInputComponent(PlayerInputComponent);

	PlayerInputComponent->BindAxis(TEXT("Throttle"), Control, &UVehicleControl::OnThrottle);
	PlayerInputComponent->BindAxis(TEXT("Brake"), Control, &UVehicleControl::OnBrake);
	PlayerInputComponent->BindAxis(TEXT("Steer"), Control, &UVehicleControl::OnSteer);
	//PlayerInputComponent->BindAxis(TEXT("Handbrake"), &Control, &FVehicleControl::OnHandbrake);

	PlayerInputComponent->BindAction(TEXT("SteerLeft"), EInputEvent::IE_Pressed, Control, &UVehicleControl::OnSteerLeftPressed);
	PlayerInputComponent->BindAction(TEXT("SteerLeft"), EInputEvent::IE_Released, Control, &UVehicleControl::OnSteerLeftReleased);
	PlayerInputComponent->BindAction(TEXT("SteerRight"), EInputEvent::IE_Pressed, Control, &UVehicleControl::OnSteerRightPressed);
	PlayerInputComponent->BindAction(TEXT("SteerRight"), EInputEvent::IE_Released, Control, &UVehicleControl::OnSteerRightReleased);

	//PlayerInputComponent->BindAction(TEXT("GearUp"), EInputEvent::IE_Pressed, &Control, &FVehicleControl::OnGearUp);
	//PlayerInputComponent->BindAction(TEXT("GearDown"), EInputEvent::IE_Pressed, &Control, &FVehicleControl::OnGearDown);
}

UVehicleControl::UVehicleControl()
	: Throttle(0.0f)
	, Brake(0.0f)
	, Steer(0.0f)
	, Gear(0)
	, bSteerByButton(false)
	, SteerAxisInput(0.0f)
	, bSteerLeft(false)
	, bSteerRight(false)
{

}

void UVehicleControl::Tick(float DeltaSeconds)
{
	TickButtonSteering(DeltaSeconds);

	if (GEngine)
	{
		GEngine->AddOnScreenDebugMessage((uint64)this, 0, FColor::White,
			FString::Printf(TEXT("Gear(%d), Steer(%04.2f) Throttle(%04.2f) Brake(%04.2f)"),
				Gear, Steer, Throttle, Brake));
	}
}

void UVehicleControl::SetThrottle(float InThrottle)
{
	ensureMsgf(InThrottle >= 0.0f && InThrottle <= 1.0f, TEXT("Invalid throttle for vehicle: %f"), InThrottle);

	Throttle = FMath::Clamp(InThrottle, 0.0f, 1.0f);
}

void UVehicleControl::SetBrake(float InBrake)
{
	ensureMsgf(InBrake >= 0.0f && InBrake <= 1.0f, TEXT("Invalid brake for vehicle: %f"), InBrake);

	Brake = FMath::Clamp(InBrake, 0.0f, 1.0f);
}

void UVehicleControl::SetSteer(float InSteer)
{
	ensureMsgf(InSteer >= -1.0f && InSteer <= 1.0f, TEXT("Invalid steer for vehicle: %f"), InSteer);

	Steer = FMath::Clamp(InSteer, -1.0f, 1.0f);
}

void UVehicleControl::OnThrottle(float InThrottle)
{
	if (InThrottle < 0)
	{
		// This can happen with pad(using same axis with brake)
		return;
	}

	SetThrottle(InThrottle);
}

void UVehicleControl::OnBrake(float InBrake)
{
	if (InBrake < 0)
	{
		// This can happen with pad(using same axis with throttle)
		return;
	}

	SetBrake(InBrake);
}

void UVehicleControl::OnSteer(float InSteer)
{
	SteerAxisInput = InSteer;

	if (bSteerByButton && InSteer != 0.0f)
	{
		bSteerByButton = false;
	}
}

void UVehicleControl::OnSteerLeftPressed()
{
	bSteerLeft = true;
	bSteerByButton = true;
}

void UVehicleControl::OnSteerLeftReleased()
{
	bSteerLeft = false;
	bSteerByButton = true;
}

void UVehicleControl::OnSteerRightPressed()
{
	bSteerRight = true;
	bSteerByButton = true;
}

void UVehicleControl::OnSteerRightReleased()
{
	bSteerRight = false;
	bSteerByButton = true;
}

void UVehicleControl::TickButtonSteering(float DeltaSeconds)
{
	if (!bSteerByButton)
	{
		return;
	}

	const float SteeringSpeed = bSteerLeft ? 3.0f : -3.0f;

	float NewSteering = GetSteer();

	if (bSteerLeft || bSteerRight)
	{
		NewSteering -= SteeringSpeed * DeltaSeconds;
	}
	else
	{
		float ReturnSteering = FMath::Min(FMath::Abs(NewSteering), 3.0f * DeltaSeconds);
		NewSteering += FMath::Sign(NewSteering) * -1 * ReturnSteering;
	}

	SetSteer(FMath::Clamp(NewSteering, -1.0f, 1.0f));
}
