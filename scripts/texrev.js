/// Look for .tex files in this repo. If the file contains the line %DRAFT_ID_MARKER%
/// replace the next non-empty line with \def\revid{md5}, where md5 is the md5 hash of the file
/// excluding the revid line
//

const { promises: fs } = require('fs')
const path = require('path')
const crypto = require('crypto')

const workdir = path.join(__dirname, "..")

const MARKER = '%DRAFT_ID_MARKER%'

hash_dir(workdir).catch(console.log)

async function hash_dir(dirpath) {
    const files = await fs.readdir(dirpath)
    for (name of files) {
        if (["target", ".git"].indexOf(name) !== -1) {
            continue
        }
        const filepath = path.join(dirpath, name)
        const stat = await fs.stat(filepath)
        if (stat.isDirectory()) {
            await hash_dir(filepath)
        } else if (name.endsWith('.tex')) {
            await update_hash(filepath)
        }
    }
}

async function update_hash(texpath) {
    const contents = (await fs.readFile(texpath)).toString()
    const lines = contents.split('\n')
    const marker_index = lines.indexOf(MARKER)
    if (marker_index === -1) {
        console.log('not found')
        return
    }
    const hash = crypto.createHash('MD5')
    hash.update(contents)
    const hex_digest = hash.digest('hex').slice(0, 8)
    lines.splice(marker_index + 1, 1, `\\def\\revid\{Draft revision ${hex_digest}\}`)
    await fs.writeFile(texpath, lines.join('\n'))
}

