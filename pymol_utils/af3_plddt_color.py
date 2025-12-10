# For open source Pymol installations, please refer to https://github.com/schrodinger/pymol-open-source

from pymol import cmd

# Define custom pLDDT colors (RGB in 0–1 range)
cmd.set_color("plddt_vh", [0.051, 0.341, 0.827])  # Very high (>90)
cmd.set_color("plddt_h",  [0.416, 0.796, 0.945])  # High (90 >= x > 70)
cmd.set_color("plddt_l",  [0.996, 0.851, 0.212])  # Low  (70 >= x > 50)
cmd.set_color("plddt_vl", [0.992, 0.490, 0.302])  # Very low (<=50)

def af3_color_plddt(selection="all"):
    """
    Color AlphaFold(3) structures by pLDDT (stored in B-factor).
    Usage: af3_color_plddt sele
    """

    # Use 'between a, b' (inclusive) to avoid '<=' parser issues
    cmd.color("plddt_vh", f"({selection}) and b>90")
    cmd.color("plddt_h",  f"({selection}) and (b>70 and not b>90)")  # ≈ 70<b<=90
    cmd.color("plddt_l",  f"({selection}) and (b>50 and not b>70)")  # ≈ 50<b<=70
    cmd.color("plddt_vl", f"({selection}) and not b>50")             # ≈ b<=50

    '''
    cmd.color("plddt_vh", f"({selection}) and b > 90")
    cmd.color("plddt_h",  f"({selection}) and b <= 90 and b > 70")
    cmd.color("plddt_l",  f"({selection}) and b <= 70 and b > 50")
    cmd.color("plddt_vl", f"({selection}) and b <= 50")
    '''

# register PyMOL command
cmd.extend("af3_color_plddt", af3_color_plddt)
cmd.auto_arg[0]["af3_color_plddt"] = [cmd.object_sc, "object", ""]