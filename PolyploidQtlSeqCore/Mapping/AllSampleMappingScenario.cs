using Kurukuru;
using PolyploidQtlSeqCore.IO;
using PolyploidQtlSeqCore.Share;

namespace PolyploidQtlSeqCore.Mapping
{
    /// <summary>
    /// 全サンプルのMappingを行うシナリオ
    /// </summary>
    internal class AllSampleMappingScenario
    {
        private readonly SampleMappingService _sampleMappingService;

        /// <summary>
        /// 全サンプルのMappingシナリオを作成する。
        /// </summary>
        /// <param name="refSeq">リファレンスシークエンス</param>
        /// <param name="thread">スレッド数</param>
        public AllSampleMappingScenario(ReferenceSequence refSeq, ThreadNumber thread)
        {
            _sampleMappingService = new SampleMappingService(refSeq, thread);
        }

        /// <summary>
        /// 全サンプルのMappingを行う。
        /// </summary>
        /// <param name="option">Mappingオプション</param>
        /// <returns>全サンプルのBAMファイル</returns>
        public async ValueTask<AllSampleBamFiles> MappingAsync(MappingSampleSettings option)
        {
            var p1BamFile = await MappingAsync(option.Parent1Directory);
            var p2BamFile = await MappingAsync(option.Parent2Directory);
            var bulk1BamFile = await MappingAsync(option.Bulk1Directory);
            var bulk2BamFile = await MappingAsync(option.Bulk2Directory);

            Log.Clear();

            return new AllSampleBamFiles(p1BamFile, p2BamFile, bulk1BamFile, bulk2BamFile);
        }

        private async ValueTask<BamFile> MappingAsync(ISampleDirectory sampleDirectory)
        {
            return await Spinner.StartAsync($"{sampleDirectory.SampleName} mapping ...", async spinner =>
            {
                var skip = false;

                try
                {
                    Log.Clear();

                    if (sampleDirectory.HasBamFile())
                    {
                        skip = true;
                        spinner.Succeed($"{sampleDirectory.SampleName} mapping skip");
                        return sampleDirectory.ToBamFile();
                    }

                    var bamFile = await _sampleMappingService.MappingAsync(sampleDirectory);
                    spinner.Succeed($"{sampleDirectory.SampleName} mapping completed");

                    return bamFile;
                }
                catch
                {
                    spinner.Fail($"{sampleDirectory.SampleName} mapping error");

                    throw;
                }
                finally
                {
                    if (!skip)
                    {
                        var logFilePath = Path.Combine(sampleDirectory.Path, $"{sampleDirectory.SampleName} Mapping.log.txt");
                        Log.Save(logFilePath);
                    }
                }
            });
        }
    }
}
