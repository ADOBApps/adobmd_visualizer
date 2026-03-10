-- ADOBMD Fast Parser for LuaJIT
-- Optimized for speed with minimal allocations
-- This file is part of Quantum Analysis Helper.
-- Quantum Analysis Helper is free software: you can redistribute it and/or modify
-- it under the terms of the GNU General Public License as published by
-- the Free Software Foundation, either version 3 of the License, or
-- (at your option) any later version.

-- Copyright (C) [2026] Acxel David Orozco Baldomero

local ffi = require("ffi")

-- C types for fast number handling
ffi.cdef[[
    typedef struct { double x, y, z; } vec3_t;
    typedef struct { int id, type, mol; double charge; int qm; } atom_t;
    typedef struct { int id, type, a1, a2; } bond_t;
]]

local M = {}

-- Fast string splitting (no regex, pure Lua)
local function split_line(line)
    local parts = {}
    local idx = 1
    for token in string.gmatch(line, "%S+") do
        parts[idx] = token
        idx = idx + 1
    end
    return parts
end

-- Parse header in ONE pass
function M.parse_header(content)
    local header = {
        natoms = 0,
        nbonds = 0,
        nangles = 0,
        ndihedrals = 0,
        nimpropers = 0,
        natom_types = 0,
        nbond_types = 0,
        box_lo = {0, 0, 0},
        box_hi = {0, 0, 0},
        has_box = false
    }
    
    for line in string.gmatch(content, "[^\n]+") do
        if string.find(line, "atoms$") then
            header.natoms = tonumber(string.match(line, "^(%d+)"))
        elseif string.find(line, "bonds$") then
            header.nbonds = tonumber(string.match(line, "^(%d+)"))
        elseif string.find(line, "atom types$") then
            header.natom_types = tonumber(string.match(line, "^(%d+)"))
        elseif string.find(line, "xlo xhi") then
            local lo, hi = string.match(line, "([%d%.%-]+)%s+([%d%.%-]+)")
            header.box_lo[1] = tonumber(lo)
            header.box_hi[1] = tonumber(hi)
            header.has_box = true
        elseif string.find(line, "ylo yhi") then
            local lo, hi = string.match(line, "([%d%.%-]+)%s+([%d%.%-]+)")
            header.box_lo[2] = tonumber(lo)
            header.box_hi[2] = tonumber(hi)
        elseif string.find(line, "zlo zhi") then
            local lo, hi = string.match(line, "([%d%.%-]+)%s+([%d%.%-]+)")
            header.box_lo[3] = tonumber(lo)
            header.box_hi[3] = tonumber(hi)
        end
    end
    
    return header
end

-- Parse masses section FAST
function M.parse_masses(content)
    local masses = {}
    local type_elements = {}
    
    local in_section = false
    for line in string.gmatch(content, "[^\n]+") do
        if line == "Masses" then
            in_section = true
        elseif in_section then
            if line == "" or line:sub(1,1) == '#' then
                -- skip comments
            elseif line:match("^%d") then
                -- Parse mass line
                local type_id = tonumber(line:match("^(%d+)"))
                local mass = tonumber(line:match("%s+(%d+%.?%d*)"))
                masses[type_id] = mass
                
                -- Extract element from comment
                local comment = line:match("#%s*(.*)")
                if comment then
                    local elem = comment:match("([A-Z][a-z]?)")
                    if elem then
                        type_elements[type_id] = elem
                    end
                end
            else
                break
            end
        end
    end
    
    return masses, type_elements
end

-- Parse atoms section - THE BOTTLENECK!
function M.parse_atoms(content, type_elements)
    local atoms = {}
    local qm_indices = {}
    local count = 0
    
    local in_section = false
    for line in string.gmatch(content, "[^\n]+") do
        if line == "Atoms" then
            in_section = true
        elseif in_section then
            if line == "" or line:sub(1,1) == '#' then
                -- skip
            elseif line:match("^%d") then
                -- FAST parsing: split once, convert once
                local parts = {}
                for token in string.gmatch(line, "%S+") do
                    parts[#parts+1] = token
                end
                
                if #parts >= 7 then
                    count = count + 1
                    local id = tonumber(parts[1])
                    local mol = tonumber(parts[2])
                    local type_id = tonumber(parts[3])
                    local charge = tonumber(parts[4])
                    local x = tonumber(parts[5])
                    local y = tonumber(parts[6])
                    local z = tonumber(parts[7])
                    
                    local is_qm = false
                    if #parts >= 8 then
                        is_qm = (parts[8] == "QM" or parts[8] == "Q")
                    end
                    
                    local element = type_elements[type_id] or "X"
                    
                    -- Store as table (fast in LuaJIT)
                    atoms[count] = {
                        id = id,
                        type = type_id,
                        mol = mol,
                        element = element,
                        x = x,
                        y = y,
                        z = z,
                        charge = charge,
                        qm = is_qm and 1 or 0
                    }
                    
                    if is_qm then
                        qm_indices[#qm_indices+1] = id
                    end
                end
            else
                break
            end
        end
    end
    
    return atoms, qm_indices
end

-- Parse bonds section
function M.parse_bonds(content)
    local bonds = {}
    local count = 0
    
    local in_section = false
    for line in string.gmatch(content, "[^\n]+") do
        if line == "Bonds" then
            in_section = true
        elseif in_section then
            if line == "" or line:sub(1,1) == '#' then
                -- skip
            elseif line:match("^%d") then
                local parts = {}
                for token in string.gmatch(line, "%S+") do
                    parts[#parts+1] = token
                end
                
                if #parts >= 4 then
                    count = count + 1
                    bonds[count] = {
                        id = tonumber(parts[1]),
                        type = tonumber(parts[2]),
                        a1 = tonumber(parts[3]),
                        a2 = tonumber(parts[4]),
                        order = tonumber(parts[2])  -- use type as order
                    }
                end
            else
                break
            end
        end
    end
    
    return bonds
end

-- MAIN FUNCTION: Parse entire file in ONE PASS
function M.parse_file(content)
    -- Single pass through file
    local header = M.parse_header(content)
    local masses, type_elements = M.parse_masses(content)
    local atoms, qm_indices = M.parse_atoms(content, type_elements)
    local bonds = M.parse_bonds(content)
    
    -- Build result as a single table
    local result = {
        header = {
            natoms = header.natoms,
            nbonds = header.nbonds,
            natom_types = header.natom_types,
            has_box = header.has_box,
            box_lo = {header.box_lo[1], header.box_lo[2], header.box_lo[3]},
            box_hi = {header.box_hi[1], header.box_hi[2], header.box_hi[3]}
        },
        type_elements = type_elements,
        atoms = atoms,
        bonds = bonds,
        qm_indices = qm_indices,
        has_qm = (#qm_indices > 0)
    }
    
    -- Convert to JSON string
    return M.to_json(result)
end

-- Simple JSON encoder
function M.to_json(obj)
    local function encode(val)
        local t = type(val)
        
        if t == "string" then
            return '"' .. val:gsub('["\\]', '\\%0') .. '"'
        elseif t == "number" then
            return tostring(val)
        elseif t == "boolean" then
            return val and "true" or "false"
        elseif t == "table" then
            -- Check if array
            local is_array = true
            local max_idx = 0
            for k in pairs(val) do
                if type(k) ~= "number" or k < 1 then
                    is_array = false
                    break
                end
                max_idx = math.max(max_idx, k)
            end
            
            if is_array and max_idx > 0 then
                -- Array
                local parts = {}
                for i = 1, max_idx do
                    parts[i] = encode(val[i]) or "null"
                end
                return "[" .. table.concat(parts, ",") .. "]"
            else
                -- Object
                local parts = {}
                for k, v in pairs(val) do
                    parts[#parts+1] = encode(tostring(k)) .. ":" .. encode(v)
                end
                return "{" .. table.concat(parts, ",") .. "}"
            end
        else
            return "null"
        end
    end
    
    return encode(obj)
end

-- CSV export function (fast)
function M.export_to_csv(atoms, bonds)
    local csv = {}
    csv[1] = "id,element,x,y,z,charge,qm"
    
    for i, atom in ipairs(atoms) do
        csv[i+1] = string.format("%d,%s,%.6f,%.6f,%.6f,%.3f,%d",
            atom.id, atom.element, atom.x, atom.y, atom.z, atom.charge, atom.qm)
    end
    
    if bonds and #bonds > 0 then
        csv[#csv+1] = ""
        csv[#csv+1] = "# bonds"
        csv[#csv+1] = "id,type,atom1,atom2"
        for i, bond in ipairs(bonds) do
            csv[#csv+1] = string.format("%d,%d,%d,%d",
                bond.id, bond.type, bond.a1, bond.a2)
        end
    end
    
    return table.concat(csv, "\n")
end

-- Return module
return M