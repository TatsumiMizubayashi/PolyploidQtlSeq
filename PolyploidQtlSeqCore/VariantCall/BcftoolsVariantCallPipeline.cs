using Cysharp.Diagnostics;
using PolyploidQtlSeqCore.IO;
using PolyploidQtlSeqCore.Mapping;
using PolyploidQtlSeqCore.QtlAnalysis;
using PolyploidQtlSeqCore.QtlAnalysis.Chr;
using static Zx.Env;

namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// Bcftools変異検出パイプライン
    /// </summary>
    internal class BcftoolsVariantCallPipeline
    {
        private readonly BcfToolsVariantCallSettings _settings;
        private const int MAX_DEPTH = 10000;

        /// <summary>
        /// Bcftools変異検出パイプラインを作成する。
        /// </summary>
        /// <param name="settings">設定</param>
        /// <param name="refSeq">リファレンスシークエンス</param>
        public BcftoolsVariantCallPipeline(BcfToolsVariantCallSettings settings)
        {
            _settings = settings;
        }

        /// <summary>
        /// 指定染色体の変異検出を行う。
        /// </summary>
        /// <param name="bamFiles">BAMファイル</param>
        /// <param name="targetChr">染色体</param>
        /// <returns>VCFファイル</returns>
        public async ValueTask<OneChromosomeVcfFile> CallAsync(AllSampleBamFiles bamFiles, Chromosome targetChr)
        {
            var outputVcfPath = _settings.OutputDirectory.CreateFilePath(targetChr.Name + ".vcf.gz");

            var command = $"bcftools mpileup -a AD,ADF,ADR -B "
                + $"-q {_settings.MinmumMappingQuality.Value} -Q {_settings.MinmumBaseQuality.Value} -C {_settings.AdjustMappingQuality.Value} "
                + $"-f {_settings.ReferenceSequence.Path} -d {MAX_DEPTH} -r {targetChr.ToRegion()} -O u "
                + $"{bamFiles.Parent1BamFile.Path} {bamFiles.Parent2BamFile.Path} {bamFiles.Bulk1BamFile.Path} {bamFiles.Bulk2BamFile.Path} "
                + "| bcftools call -v -m -a GQ,GP -O u "
                + $"| bcftools filter -i 'INFO/MQ>={_settings.MinmumMappingQuality.Value}' -O u "
                + $"| bcftools norm -f {_settings.ReferenceSequence.Path} -O u "
                + $"| bcftools sort -O z -o {outputVcfPath}";
            CommandLog.Add(command);

            try
            {
                verbose = false;
                var (_, stdErrors) = await processl2(command);
                if (stdErrors.Length != 0)
                {
                    Log.Add($"{targetChr.Name} Variant Call");
                    Log.AddRange(stdErrors);
                    Log.AddSeparator();
                }

                var vcfFile = new OneChromosomeVcfFile(outputVcfPath, targetChr);
                await vcfFile.CreateIndexFile();

                return vcfFile;
            }
            catch (ProcessErrorException ex)
            {
                Log.Add(ex);

                throw;
            }
        }
    }
}
