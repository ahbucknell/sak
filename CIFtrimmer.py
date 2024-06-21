from gemmi import cif
import more_itertools as mit

doc = cif.read_file('test.cif')
block = doc.sole_block()
atom_loop = block.find_loop('_atom_site.id')
amino_loop = block.find_loop('_atom_site.label_seq_id')
plddt_loop = block.find_loop('_atom_site.B_iso_or_equiv')

def avg_pLDDT(atoms, aminos, plddt):
    total = []
    locs = []
    avg = {}
    maps = {}
    for i in range(len(aminos)):
        now = int(aminos[i])
        if i + 1 < len(aminos):
            after = int(aminos[i + 1])
            score = float(plddt[i])
            total.append(score)
            locs.append(int(atoms[i]))
            if after != now:
                avg[now] = sum(total) / len(total)
                maps[now] = locs
                total = []
                locs = []
        else:
            end = True
            score = float(plddt[i])
            total.append(score)
            locs.append(int(atoms[i]))
            avg[now] = sum(total) / len(total)
            maps[now] = locs
    return avg, maps

def trimmer(avg_pLDDT, threshold):
    for i in range(1, len(avg_pLDDT)+1):
        if float(avg_pLDDT[i]) < int(threshold):
            avg_pLDDT.pop(i)
    return avg_pLDDT

def find_ranges(iterable):
    for group in mit.consecutive_groups(iterable):
        group = list(group)
        if len(group) > 1:
            yield group[0], group[-1]
            
def trim_ranges(ranges, threshold):
    trimranges = []
    finalAAs = []
    for i in range(0, len(ranges)):
        delta = (ranges[i][1]-ranges[i][0])+1
        if delta >= threshold:
            trimranges.append(ranges[i])
    for j in range(0,len(trimranges)):
        for k in range(trimranges[j][0],trimranges[j][1]+1):
            finalAAs.append(k)
    return finalAAs

def amino_atom_conv(aminos, atom_map):
    final_atoms = []
    for i in range(0, len(aminos)):
        final_atoms.append(atom_map[aminos[i]])
    return list(mit.flatten(final_atoms))    

avg, keymap = avg_pLDDT(atom_loop, amino_loop, plddt_loop)
trimmed = trimmer(avg, 70)
consecutive = list(find_ranges(trimmed))
final_aminos = trim_ranges(consecutive, 50)
final_atoms = amino_atom_conv(final_aminos, keymap)

newdoc = cif.Document()
newblock = newdoc.add_new_block(str(block.name))
newblock.set_pair('_entry.id', newblock.name)
newblock.set_pair('_audit_confirm.dict_location', block.find_pair('_audit_conform.dict_location')[1])
newblock.set_pair('_audit_confirm.dict_name', block.find_pair('_audit_conform.dict_name')[1])
newblock.set_pair('_audit_confirm.dict_version', block.find_pair('_audit_conform.dict_version')[1])

entity = newblock.init_loop('_entity.',['id','pdbx_description','type'])
for a in range(0, len(block.find_loop("_entity.id"))):
    entity.add_row([
        block.find_loop("_entity.id")[a],
        block.find_loop("_entity.pdbx_description")[a],
        block.find_loop("_entity.type")[a]
        ])

entity_poly = newblock.init_loop('_entity_poly.',['entity_id', 'pdbx_strand_id', 'type'])
for b in range(0, len(block.find_loop('_entity_poly.entity_id'))):
    entity_poly.add_row([
        block.find_loop("_entity_poly.entity_id")[b],
        block.find_loop("_entity_poly.pdbx_strand_id")[b],
        block.find_loop("_entity_poly.type")[b]
        ])
    
entity_poly_seq = newblock.init_loop('_entity_poly_seq.',['entity_id', 'hetero', 'mon_id','num'])
for c in range(0, len(block.find_loop('_entity_poly_seq.entity_id'))):
    entity_poly_seq.add_row([
        block.find_loop("_entity_poly_seq.entity_id")[c],
        block.find_loop("_entity_poly_seq.hetero")[c],
        block.find_loop("_entity_poly_seq.mon_id")[c],
        block.find_loop("_entity_poly_seq.num")[c]])

group_PDB = []
atomid = []
type_symbol = []
label_atom_id = []
label_alt_id = []
label_comp_id = []
label_asym_id = []
label_entity_id = []
label_seq_id = []
pdbx_PDB_ins_code = []
Cartn_x = []
Cartn_y = []
Cartn_z = []
occupancy = []
B_iso_or_equiv = []
auth_seq_id = []
auth_asym_id = []
pdbx_PDB_model_num = []

for x in range(0, len(final_atoms)):
    group_PDB.append(str(block.find_loop("_atom_site.group_PDB")[final_atoms[x]-1]))
    atomid.append(str(block.find_loop("_atom_site.id")[final_atoms[x]-1]))
    type_symbol.append(str(block.find_loop("_atom_site.type_symbol")[final_atoms[x]-1]))
    label_atom_id.append(str(block.find_loop("_atom_site.label_atom_id")[final_atoms[x]-1]))
    label_alt_id.append(str(block.find_loop("_atom_site.label_alt_id")[final_atoms[x]-1]))
    label_comp_id.append(str(block.find_loop("_atom_site.label_comp_id")[final_atoms[x]-1]))
    label_asym_id.append(str(block.find_loop("_atom_site.label_asym_id")[final_atoms[x]-1]))
    label_entity_id.append(str(block.find_loop("_atom_site.label_entity_id")[final_atoms[x]-1]))
    label_seq_id.append(str(block.find_loop("_atom_site.label_seq_id")[final_atoms[x]-1]))
    pdbx_PDB_ins_code.append(str(block.find_loop("_atom_site.pdbx_PDB_ins_code")[final_atoms[x]-1]))
    Cartn_x.append(str(block.find_loop("_atom_site.Cartn_x")[final_atoms[x]-1]))
    Cartn_y.append(str(block.find_loop("_atom_site.Cartn_y")[final_atoms[x]-1]))
    Cartn_z.append(str(block.find_loop("_atom_site.Cartn_z")[final_atoms[x]-1]))
    occupancy.append(str(block.find_loop("_atom_site.occupancy")[final_atoms[x]-1]))
    B_iso_or_equiv.append(str(block.find_loop("_atom_site.B_iso_or_equiv")[final_atoms[x]-1]))
    auth_seq_id.append(str(block.find_loop("_atom_site.auth_seq_id")[final_atoms[x]-1]))
    auth_asym_id.append(str(block.find_loop("_atom_site.auth_asym_id")[final_atoms[x]-1]))
    pdbx_PDB_model_num.append(str(block.find_loop("_atom_site.pdbx_PDB_model_num")[final_atoms[x]-1]))

atom_site = newblock.init_loop('_atom_site.',
                          ['group_PDB', 'id','type_symbol','label_atom_id','label_alt_id','label_comp_id',
                           'label_asym_id','label_entity_id','label_seq_id','pdbx_PDB_ins_code','Cartn_x',
                           'Cartn_y','Cartn_z','occupancy','B_iso_or_equiv','auth_seq_id','auth_asym_id',
                           'pdbx_PDB_model_num'])
for y in range(0, len(group_PDB)):
    atom_site.add_row(cif.quote_list([
        group_PDB[y],
        atomid[y],
        type_symbol[y],
        label_atom_id[y],
        label_alt_id[y],
        label_comp_id[y],
        label_asym_id[y],
        label_entity_id[y],
        label_seq_id[y],
        pdbx_PDB_ins_code[y],
        Cartn_x[y],
        Cartn_y[y],
        Cartn_z[y],
        occupancy[y],
        B_iso_or_equiv[y],
        auth_seq_id[y],
        auth_asym_id[y],
        pdbx_PDB_model_num[y]]))
newdoc.write_file('test-trimmed.cif')
print("Done!")
