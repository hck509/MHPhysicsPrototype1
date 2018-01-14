// Fill out your copyright notice in the Description page of Project Settings.

#include "VehicleActor.h"


// Sets default values
AVehicleActor::AVehicleActor()
{
 	// Set this pawn to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;

	MeshComponent = CreateDefaultSubobject<UMHStaticMeshComponent>(TEXT("Mesh"));
}

// Called when the game starts or when spawned
void AVehicleActor::BeginPlay()
{
	Super::BeginPlay();
	
}

// Called every frame
void AVehicleActor::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

}

// Called to bind functionality to input
void AVehicleActor::SetupPlayerInputComponent(UInputComponent* PlayerInputComponent)
{
	Super::SetupPlayerInputComponent(PlayerInputComponent);

	PlayerInputComponent->BindAxis(TEXT("Throttle"), &Control, &FVehicleControl::OnThrottle);
	PlayerInputComponent->BindAxis(TEXT("Brake"), &Control, &FVehicleControl::OnBrake);
	PlayerInputComponent->BindAxis(TEXT("Steer"), &Control, &FVehicleControl::OnSteer);
	//PlayerInputComponent->BindAxis(TEXT("Handbrake"), &Control, &FVehicleControl::OnHandbrake);

	PlayerInputComponent->BindAction(TEXT("SteerLeft"), EInputEvent::IE_Pressed, &Control, &FVehicleControl::OnSteerLeftPressed);
	PlayerInputComponent->BindAction(TEXT("SteerLeft"), EInputEvent::IE_Released, &Control, &FVehicleControl::OnSteerLeftReleased);
	PlayerInputComponent->BindAction(TEXT("SteerRight"), EInputEvent::IE_Pressed, &Control, &FVehicleControl::OnSteerRightPressed);
	PlayerInputComponent->BindAction(TEXT("SteerRight"), EInputEvent::IE_Released, &Control, &FVehicleControl::OnSteerRightReleased);

	//PlayerInputComponent->BindAction(TEXT("GearUp"), EInputEvent::IE_Pressed, &Control, &FVehicleControl::OnGearUp);
	//PlayerInputComponent->BindAction(TEXT("GearDown"), EInputEvent::IE_Pressed, &Control, &FVehicleControl::OnGearDown);
}

FVehicleControl::FVehicleControl()
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

void FVehicleControl::Tick(float DeltaSeconds)
{
	TickButtonSteering(DeltaSeconds);
}

void FVehicleControl::SetThrottle(float InThrottle)
{
	ensureMsgf(InThrottle >= 0.0f && InThrottle <= 1.0f, TEXT("Invalid throttle for vehicle: %f"), InThrottle);

	Throttle = FMath::Clamp(InThrottle, 0.0f, 1.0f);
}

void FVehicleControl::SetBrake(float InBrake)
{
	ensureMsgf(InBrake >= 0.0f && InBrake <= 1.0f, TEXT("Invalid brake for vehicle: %f"), InBrake);

	Brake = FMath::Clamp(InBrake, 0.0f, 1.0f);
}

void FVehicleControl::SetSteer(float InSteer)
{
	ensureMsgf(InSteer >= -1.0f && InSteer <= 1.0f, TEXT("Invalid steer for vehicle: %f"), InSteer);

	Steer = FMath::Clamp(InSteer, -1.0f, 1.0f);
}

void FVehicleControl::OnThrottle(float InThrottle)
{
	if (InThrottle < 0)
	{
		// This can happen with pad(using same axis with brake)
		return;
	}

	SetThrottle(InThrottle);
}

void FVehicleControl::OnBrake(float InBrake)
{
	if (InBrake < 0)
	{
		// This can happen with pad(using same axis with throttle)
		return;
	}

	SetBrake(InBrake);
}

void FVehicleControl::OnSteer(float InSteer)
{
	SteerAxisInput = InSteer;

	if (bSteerByButton && InSteer != 0.0f)
	{
		bSteerByButton = false;
	}
}

void FVehicleControl::OnSteerLeftPressed()
{
	bSteerLeft = true;
	bSteerByButton = true;
}

void FVehicleControl::OnSteerLeftReleased()
{
	bSteerLeft = false;
	bSteerByButton = true;
}

void FVehicleControl::OnSteerRightPressed()
{
	bSteerRight = true;
	bSteerByButton = true;
}

void FVehicleControl::OnSteerRightReleased()
{
	bSteerRight = false;
	bSteerByButton = true;
}

void FVehicleControl::TickButtonSteering(float DeltaSeconds)
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
