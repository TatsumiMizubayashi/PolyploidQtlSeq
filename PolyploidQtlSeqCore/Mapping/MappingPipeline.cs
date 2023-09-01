using Cysharp.Diagnostics;
using PolyploidQtlSeqCore.IO;
using static Zx.Env;

namespace PolyploidQtlSeqCore.Mapping
{
    /// <summary>
    /// Mappingパイプライン
    /// </summary>
    internal class MappingPipeline
    {
        private readonly MappingSettings _settings;

        /// <summary>
        /// Mappingパイプラインを作成する。
        /// </summary>
        /// <param name="settings">Mapping設定</param>
        public MappingPipeline(MappingSettings settings)
        {
            _settings = settings;
        }

        /// <summary>
        /// Mappingを行う。
        /// </summary>
        /// <param name="sampleName">サンプル名</param>
        /// <param name="fastqFilePair">Fastqファイルペア</param>
        /// <returns>BAMファイル</returns>
        public async ValueTask<BamFile> MappingAsync(string sampleName, FastqFilePair fastqFilePair)
        {
            var readGroup = CreateReadGroup(sampleName, fastqFilePair);
            var bamFilePath = CreateBamFilePath(fastqFilePair);

            var refSeq = _settings.ReferenceSequence.Path;
            var threadNum = _settings.ThreadNumber.Value;
            var command = $"bwa mem -t {threadNum} -M -R '{readGroup}' {refSeq} {fastqFilePair.Fastq1Path} {fastqFilePair.Fastq2Path} "
                + "| samtools fixmate -m - - "
                + $"| samtools sort -@ {threadNum} "
                + "| samtools markdup -r - - "
                + $"| samtools view -b -f 2 -F 2048 -o {bamFilePath}";
            CommandLog.Add(command);

            try
            {
                verbose = false; // コンソールに出力しない
                var (_, stdErrors) = await processl2(command);
                Log.AddRange(stdErrors);
            }
            catch (ProcessErrorException ex)
            {
                Log.Add(ex);

                throw;
            }

            var bamFile = new BamFile(sampleName, bamFilePath);
            await bamFile.CreateIndexFileAsync();

            return bamFile;
        }

        /// <summary>
        /// リードグループを作成する。
        /// </summary>
        /// <param name="sampleName">サンプル名</param>
        /// <param name="fastqFilePair">Fastqファイルペア</param>
        /// <returns>リードグループ</returns>
        private static string CreateReadGroup(string sampleName, FastqFilePair fastqFilePair)
        {
            return @$"@RG\tID:{fastqFilePair.BaseName}\tSM:{sampleName}\tPL:illumina\tLB:{fastqFilePair.BaseName}_library";

        }

        /// <summary>
        /// BAMファイルのPathを作成する。
        /// </summary>
        /// <param name="fastqFilePair">Fastqファイルペア</param>
        /// <returns>BAMファイルPath</returns>
        private static string CreateBamFilePath(FastqFilePair fastqFilePair)
        {
            var dirPath = Path.GetDirectoryName(fastqFilePair.Fastq1Path) ?? "";

            return Path.Combine(dirPath, fastqFilePair.BaseName + ".bam");
        }
    }
}
