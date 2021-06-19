using BVHFiles
using Test


@testset "file" begin
    g = load("Example.bvh") |>
        global_positions! |>
        remove_joint!(7) |>
        remove_joint!(13) |>
        remove_joints!("J_L_Bale", "J_R_Bale", "J_L4", "J_L3", "J_L1", "J_T12", "J_T10", "J_T9", "J_T8", "J_T6", "J_T5", "J_T4", "J_T3", "J_T2") |>
        remove_joint!("J_T1", "J_C7") |>
        remove_joint!("J_C7", "J_C6") |>
        remove_joint!("J_C6", "J_C5") |>
        remove_joints!("J_C4", "J_C3", "J_C2", "J_C1", "J_Atlas") |>
        optimize_offsets!

    @test nv(g) == 25
    @test ne(g) == 24

    dict = Dict(  "J_L_Hip" => "lThighBend", 
                    "J_R_Hip" => "rThighBend", 
                    "J_L_Knee" => "lShin", 
                    "J_R_Knee" => "rShin", 
                    "J_L_Ankle" => "lFoot", 
                    "J_R_Ankle" => "rFoot", 
                    "J_L5" => "abdomenLower", 
                    "J_L2" => "abdomenUpper", 
                    "J_T11" => "chestLower", 
                    "J_T7" => "chestUpper", 
                    "J_L_Clavicle" => "lCollar", 
                    "J_R_Clavicle" => "rCollar", 
                    "J_L_Shoulder" => "lShldrBend", 
                    "J_R_Shoulder" => "rShldrBend", 
                    "J_L_Elbow" => "lForearmBend", 
                    "J_R_Elbow" => "rForearmBend", 
                    "J_L_Hand" => "lHand", 
                    "J_R_Hand" => "rHand", 
                    "J_C5" => "neckLower")

    g |>
        add_joint!("J_L_Hip", "J_L_Knee", "lThighTwist") |>
        add_joint!("J_R_Hip", "J_R_Knee", "rThighTwist") |>
        add_joint!("J_L_Shoulder", "J_L_Elbow", "lShldrTwist") |>
        add_joint!("J_R_Shoulder", "J_R_Elbow", "rShldrTwist") |>
        add_joint!("J_L_Elbow", "J_L_Hand", "lForearmTwist") |>
        add_joint!("J_R_Elbow", "J_R_Hand", "rForearmTwist") |>
        rename!(dict)

    @test nv(g) == 31
    @test ne(g) == 30

    name!(g, 1, "ROOT hip")
    oabdomen = offset(g, 1, find(g, "abdomenLower"))
    add_joint!(g, "hip", "pelvis", [0.0, 0.0, oabdomen[3]], [find(g, "lThighBend"), find(g, "rThighBend")])
    ofoot = offset(g, find(g, "lShin"), find(g, "lFoot"))
    add_joint!(g, "lFoot", "lMetatarsals", [0.0, ofoot[3] / 70, ofoot[3] / 70])
    add_joint!(g, "rFoot", "rMetatarsals", [0.0, ofoot[3] / 70, ofoot[3] / 70])
    exclude = [find(g, "lCollar"), find(g, "rCollar"), find(g, "lThighBend"), find(g, "rThighBend")]
    T = [-1.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0]

    @test nv(g) == 34
    @test ne(g) == 33
end