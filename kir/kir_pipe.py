import subprocess
import uuid
import re
import glob
from typing import ClassVar, Iterable, Type, Any


class FileMod:
    def getID(self, name: str) -> str:
        """Get id from filename"""
        # TODO: Fix id extraction
        # input_name.template_args[0]
        print(name)
        id = re.findall(r"\.(\d+)", name)[0].strip()
        return str(id)

    def listFiles(self, name: str) -> list[str]:
        """List the file names that match the wildcard"""
        new_name_set = set()
        for possible_name in glob.glob(name.replace("{}", "*") + "*"):
            # print(r"([^\.]*)".join(map(re.escape, name.split("{}"))), possible_name)
            ids = re.findall(
                r"([^\.]*)".join(map(re.escape, name.split("{}"))), possible_name
            )
            if not ids:
                continue
            new_name = name.format(ids[0])
            if new_name not in new_name_set:
                new_name_set.add(new_name)
        return sorted(new_name_set)

    def replaceWildcard(self, name: str, new_name: str) -> str:
        """Replace wildcard character to specific name"""
        # .replace_wildcard("_merge_depth")
        return name.replace(".{}", new_name)


class Executor:
    def __init__(self, engine_type: str = "podman") -> None:
        assert engine_type in ["podman", "docker"]
        self.engine = engine_type

    def runShell(
        self, cmd: str, cwd: str | None = None
    ) -> subprocess.CompletedProcess[str]:
        """wrap os.system"""
        print(cmd)
        proc = subprocess.run(
            cmd,
            shell=True,
            cwd=cwd,
            check=True,
            universal_newlines=True,
        )
        return proc

    def runDocker(
        self,
        image: str,
        cmd: str,
        cwd: str = "",
        opts: str = "",
    ) -> subprocess.CompletedProcess[str]:
        """run docker container"""
        random_name = str(uuid.uuid4()).split("-", 1)[0]
        # bad but works
        # if getEngine() == "singularity":
        #     opts = opts.replace(" -e ", " --env ")
        # conf = engine_config[getEngine()]
        proc = self.runShell(
            f"{self.engine} run -it --rm --name {random_name}"
            f"       {opts} -v $PWD:/app -w /app/{cwd}"
            f"       {image} {cmd}"
        )
        return proc

    def checkImage(self, image: str) -> bool:
        """Build docker container"""
        try:
            self.runShell(
                f"sh -c 'if [ ! $(podman image ls {image} -q) ]; then exit 1; fi'"
            )
            return True
        except subprocess.CalledProcessError:
            return False

    def buildImage(
        self, image: str, dockerfile: str, folder: str = ".", args: dict[str, str] = {}
    ) -> subprocess.CompletedProcess[str]:
        """Build docker container"""
        txt_arg = ""
        for key, value in args.items():
            txt_arg += f" --build-arg {key}={value} "
        return self.runShell(
            f"podman build {folder} -f {dockerfile} -t {image} {txt_arg}"
        )


class KirPipe:
    """The parent class for each KIR pipeline"""

    name: ClassVar[str] = ""

    def __init__(
        self,
        threads: int = 4,
        file_adapter: Type[FileMod] = FileMod,
        executor: Type[Executor] = Executor,
    ) -> None:
        self.images: dict[str, str] = {}
        self.file_adapter = file_adapter()
        self.executor = executor()
        self.threads = threads

    def getThreads(self) -> int:
        """Get theads for one job"""
        return self.threads

    def setThreads(self, threads: int = 1) -> None:
        """Set theads for job"""
        self.threads = threads

    def runShell(
        self, cmd: str, cwd: str | None = None
    ) -> subprocess.CompletedProcess[str]:
        """Run shell command"""
        return self.executor.runShell(cmd, cwd)

    def runDocker(
        self,
        image: str,
        cmd: str,
        cwd: str = "",
        opts: str = "",
    ) -> subprocess.CompletedProcess[str]:
        """Run docker-like command"""
        image = self.images.get(image, image)
        return self.executor.runDocker(image, cmd, cwd, opts)

    def checkImage(self, image: str) -> bool:
        """Build docker container"""
        image = self.images.get(image, image)
        return self.executor.checkImage(image)

    def buildImage(
        self, image: str, dockerfile: str, folder: str = ".", args: dict[str, str] = {}
    ) -> subprocess.CompletedProcess[str]:
        """Build docker image"""
        image = self.images.get(image, image)
        return self.executor.buildImage(image, dockerfile, folder, args)

    def samtobam(self, name: str, keep: bool = False) -> None:
        """samfile -> sorted bamfile and index (This is so useful)"""
        self.runDocker(
            "samtools", f"samtools sort  -@{self.getThreads()} {name}.sam -o {name}.bam"
        )
        self.runDocker("samtools", f"samtools index -@{self.getThreads()} {name}.bam")
        if not keep:
            self.runShell(f"rm {name}.sam")

    def getID(self, name: str) -> str:
        """Extract ID from filename"""
        return self.file_adapter.getID(name)

    def listFiles(self, name: str) -> list[str]:
        """List filename when matching the name"""
        return self.file_adapter.listFiles(name)

    def replaceWildcard(self, name: str, new_name: str) -> str:
        """Replace wildcard character `{}` in name to `new_name`"""
        return self.file_adapter.replaceWildcard(name, new_name)

    def runAll(self, input_name: str) -> str:
        """Run all the script(Don't use this when building pipeline"""
        return input_name

    def escapeName(self, name: str) -> str:
        """Replace non-word character to '-'"""
        return name.replace(".", "_").replace("/", "_")
