/// Watch all .tex files for changes and rebuild pdf automatically
//

const fs = require('fs/promises')
const path = require('path')
const child_process = require('child_process')

const bin = 'tectonic'

begin_watch(path.join(__dirname, '..')).catch(console.log)

async function begin_watch(dirpath) {
    const filenames = await fs.readdir(dirpath)
    for (name of filenames) {
        const filepath = path.join(dirpath, name)
        const stat = await fs.stat(filepath)
        if (stat.isDirectory()) {
            await begin_watch(filepath)
        } else if (name.endsWith('.tex')) {
            watch_file(filepath)
        }
    }
}

async function watch_file(filepath) {
    console.log(`watching: ${filepath}`)
    for await (const event of fs.watch(filepath)) {
        console.log(`rebuilding: ${filepath}`)
        // file at filepath changed
        const proc = child_process.fork(bin, [filepath])
        proc.on('exit', () => {
            console.log(`built: ${filepath}`)
        })
            // , (err, stdout, stderr) => {
            // if (err) {
            //     console.log(err)
            //     console.log("Error building .tex file: ${filepath}")
            //     process.exit(1)
            // }
            // if (stderr) {
            //     console.log(`${bin} stderr: ${stderr}`)
            // }
        // })
    }

}
