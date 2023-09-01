using Kurukuru;
using PolyploidQtlSeqCore.IO;

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
        /// <param name="setting">Mapping設定</param>
        public AllSampleMappingScenario(MappingSettings setting)
        {
            _sampleMappingService = new SampleMappingService(setting);
        }

        /// <summary>
        /// 全サンプルのMappingを行う。
        /// </summary>
        /// <param name="sampleSetting">Mappingサンプル設定</param>
        /// <returns>全サンプルのBAMファイル</returns>
        public async ValueTask<AllSampleBamFiles> MappingAsync(MappingSampleSettings sampleSetting)
        {
            var p1BamFile = await MappingAsync(sampleSetting.Parent1Directory);
            var p2BamFile = await MappingAsync(sampleSetting.Parent2Directory);
            var bulk1BamFile = await MappingAsync(sampleSetting.Bulk1Directory);
            var bulk2BamFile = await MappingAsync(sampleSetting.Bulk2Directory);

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
