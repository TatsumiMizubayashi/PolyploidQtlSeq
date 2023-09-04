using System.ComponentModel.DataAnnotations;
using McMaster.Extensions.CommandLineUtils;
using PolyploidQtlSeqCore.Application.QtlAnalysis;
using PolyploidQtlSeqCore.QtlAnalysis.OxyGraph;
using d = PolyploidQtlSeqCore.QtlAnalysis.Distribution;
using io = PolyploidQtlSeqCore.QtlAnalysis.IO;
using og = PolyploidQtlSeqCore.QtlAnalysis.OxyGraph;
using op = PolyploidQtlSeqCore.Options;
using qa = PolyploidQtlSeqCore.QtlAnalysis;
using qsf = PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter;
using sw = PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow;

namespace PolyploidQtlSeq
{
    /// <summary>
    /// QTL-Seq解析を行うコマンド
    /// </summary>
    [Command("qtl", Description = "QTL-Seq analysis.")]
    internal sealed class QtlAnalysisCommand : CommandBase, IQtlSeqAnalysisSettingValueOld
    {
        /// <summary>
        /// QTL解析コマンドを作成する。
        /// </summary>
        public QtlAnalysisCommand()
        {
            InputVcf = "";
            OutputDir = "";
            Parent1MostAlleleRateThreshold = qsf.Parent1MostAlleleRateThreshold.DEFAULT;
            Parent2SnpIndexRange = qsf.Parent2SnpIndexRange.DEFAULT;
            MinimumDepthThreshold = qsf.MinimumDepthThreshold.DEFAULT;
            MaxBulkSnpIndexThreshold = qsf.MaxBulkSnpIndexThreshold.DEFAULT;
            Ploidy = d.Ploidy.DEFAULT;
            Parent2PlexNumber = d.Parent2PlexNumber.DEFAULT;
            Bulk1Number = d.Bulk1Number.DEFAULT;
            Bulk2Number = d.Bulk2Number.DEFAULT;
            ReplicatesNumber = d.ReplicatesNumber.DEFAULT;
            WindowSize = sw.WindowSize.DEFAULT;
            StepSize = sw.StepSize.DEFAULT;
            FigureWidth = og.FigureWidth.DEFAULT;
            FigureHeight = og.FigureHeight.DEFAULT;
            XAxisMajorStep = og.XAxisMajorStep.DEFAULT;
            DisplayAnnotationImpacts = io.DisplayAnnotationImpacts.DEFAULT;
            ThreadNumber = PolyploidQtlSeqCore.Share.ThreadNumber.DEFAULT;
            ParameterFile = "";
        }

        [Option(ShortName = qa.InputVcf.SHORT_NAME, LongName = qa.InputVcf.LONG_NAME,
            Description = qa.InputVcf.DESCRIPTION, ValueName = "")]
        public string InputVcf { get; set; }

        [Option(ShortName = qa.OutputDirectory.SHORT_NAME, LongName = qa.OutputDirectory.LONG_NAME,
            Description = qa.OutputDirectory.DESCRIPTION, ValueName = "")]
        public string OutputDir { get; set; }

        [Option(ShortName = qsf.Parent1MostAlleleRateThreshold.SHORT_NAME, LongName = qsf.Parent1MostAlleleRateThreshold.LONG_NAME,
            Description = qsf.Parent1MostAlleleRateThreshold.DESCRIPTION, ValueName = "")]
        public double Parent1MostAlleleRateThreshold { get; set; }

        [Option(ShortName = qsf.Parent2SnpIndexRange.SHORT_NAME, LongName = qsf.Parent2SnpIndexRange.LONG_NAME,
            Description = qsf.Parent2SnpIndexRange.DESCRIPTION, ValueName = "")]
        public string Parent2SnpIndexRange { get; set; }

        [Option(ShortName = qsf.MinimumDepthThreshold.SHORT_NAME, LongName = qsf.MinimumDepthThreshold.LONG_NAME,
            Description = qsf.MinimumDepthThreshold.DESCRIPTION, ValueName = "")]
        public int MinimumDepthThreshold { get; set; }

        [Option(ShortName = qsf.MaxBulkSnpIndexThreshold.SHORT_NAME, LongName = qsf.MaxBulkSnpIndexThreshold.LONG_NAME,
            Description = qsf.MaxBulkSnpIndexThreshold.DESCRIPTION, ValueName = "")]
        public double MaxBulkSnpIndexThreshold { get; set; }


        [Option(ShortName = d.Ploidy.SHORT_NAME, LongName = d.Ploidy.LONG_NAME,
            Description = d.Ploidy.DESCRIPTION, ValueName = "")]
        [Range(d.Ploidy.MINIMUM, d.Ploidy.MAXIMUM, ErrorMessage = d.Ploidy.VALIDATION_ERROR_MESSAGE)]
        public int Ploidy { get; set; }

        [Option(ShortName = d.Parent2PlexNumber.SHORT_NAME, LongName = d.Parent2PlexNumber.LONG_NAME,
            Description = d.Parent2PlexNumber.DESCRIPTION, ValueName = "")]
        [Range(d.Parent2PlexNumber.MINIMUM, d.Parent2PlexNumber.MAXIMUM,
            ErrorMessage = d.Parent2PlexNumber.VALIDATION_ERROR_MESSAGE)]
        public int Parent2PlexNumber { get; set; }

        [Option(ShortName = d.Bulk1Number.SHORT_NAME, LongName = d.Bulk1Number.LONG_NAME,
            Description = d.Bulk1Number.DESCRIPTION, ValueName = "")]
        [Range(d.Bulk1Number.MINIMUM, d.Bulk1Number.MAXIMUM, ErrorMessage = d.Bulk1Number.VALIDATION_ERROR_MESSAGE)]
        public int Bulk1Number { get; set; }

        [Option(ShortName = d.Bulk2Number.SHORT_NAME, LongName = d.Bulk2Number.LONG_NAME,
            Description = d.Bulk2Number.DESCRIPTION, ValueName = "")]
        [Range(d.Bulk2Number.MINIMUM, d.Bulk2Number.MAXIMUM, ErrorMessage = d.Bulk2Number.VALIDATION_ERROR_MESSAGE)]
        public int Bulk2Number { get; set; }

        [Option(ShortName = d.ReplicatesNumber.SHORT_NAME, LongName = d.ReplicatesNumber.LONG_NAME,
            Description = d.ReplicatesNumber.DESCRIPTION, ValueName = "")]
        [Range(d.ReplicatesNumber.MINIMUM, d.ReplicatesNumber.MAXIMUM, ErrorMessage = d.ReplicatesNumber.VALIDATION_ERROR_MESSAGE)]
        public int ReplicatesNumber { get; set; }

        [Option(ShortName = sw.WindowSize.SHORT_NAME, LongName = sw.WindowSize.LONG_NAME,
            Description = sw.WindowSize.DESCRIPTION, ValueName = "")]
        [Range(sw.WindowSize.MINIMUM, sw.WindowSize.MAXIMUM, ErrorMessage = sw.WindowSize.VALIDATION_ERROR_MESSAGE)]
        public int WindowSize { get; set; }

        [Option(ShortName = sw.StepSize.SHORT_NAME, LongName = sw.StepSize.LONG_NAME,
            Description = sw.StepSize.DESCRIPTION, ValueName = "")]
        [Range(sw.StepSize.MINIMUM, sw.StepSize.MAXIMUM, ErrorMessage = sw.StepSize.VALIDATION_ERROR_MESSAGE)]
        public int StepSize { get; set; }

        [Option(ShortName = og.FigureWidth.SHORT_NAME, LongName = og.FigureWidth.LONG_NAME,
            Description = og.FigureWidth.DESCRIPTION, ValueName = "")]
        [Range(og.FigureWidth.MINIMUM, og.FigureWidth.MAXIMUM, ErrorMessage = og.FigureWidth.VALIDATION_ERROR_MESSAGE)]
        public int FigureWidth { get; set; }

        [Option(ShortName = og.FigureHeight.SHORT_NAME, LongName = og.FigureHeight.LONG_NAME,
            Description = og.FigureHeight.DESCRIPTION, ValueName = "")]
        [Range(og.FigureHeight.MINIMUM, og.FigureHeight.MAXIMUM, ErrorMessage = og.FigureHeight.VALIDATION_ERROR_MESSAGE)]
        public int FigureHeight { get; set; }

        [Option(ShortName = og.XAxisMajorStep.SHORT_NAME, LongName = og.XAxisMajorStep.LONG_NAME, Description = og.XAxisMajorStep.DESCRIPTION, ValueName = "")]
        [Range(og.XAxisMajorStep.MINIMUM, og.XAxisMajorStep.MAXIMUM, ErrorMessage = og.XAxisMajorStep.VALIDATION_ERROR_MESSAGE)]
        public int XAxisMajorStep { get; set; }


        [Option(ShortName = io.DisplayAnnotationImpacts.SHORT_NAME, LongName = io.DisplayAnnotationImpacts.LONG_NAME,
            Description = io.DisplayAnnotationImpacts.DESCRIPTION, ValueName = "")]
        public string DisplayAnnotationImpacts { get; set; }

        [Option(ShortName = PolyploidQtlSeqCore.Share.ThreadNumber.SHORT_NAME, LongName = PolyploidQtlSeqCore.Share.ThreadNumber.LONG_NAME,
            Description = PolyploidQtlSeqCore.Share.ThreadNumber.DESCRIPTION, ValueName = "")]
        [Range(PolyploidQtlSeqCore.Share.ThreadNumber.MINIMUM, PolyploidQtlSeqCore.Share.ThreadNumber.MAXIMUM, ErrorMessage = PolyploidQtlSeqCore.Share.ThreadNumber.VALIDATION_ERROR_MESSAGE)]
        public int ThreadNumber { get; set; }

        [Option(ShortName = op.ParameterFileParser.SHORT_NAME, LongName = op.ParameterFileParser.LONG_NAME,
            Description = op.ParameterFileParser.DESCRIPTION, ValueName = "")]
        public string ParameterFile { get; set; }


        public override async Task<int> OnExecuteAsync(CommandLineApplication app)
        {
            int code;
            try
            {
                code = await RunAsync(app.Options);
            }
            catch(Exception ex)
            {
                Console.Error.WriteLine(ex.Message);
                code = 1;
            }

            return code;
        }

        // OnExecuteが非同期なのでTaskを使用して処理する
        private async ValueTask<int> RunAsync(IReadOnlyCollection<CommandOption> options)
        {
            return await Task.Run(() =>
            {
                var analyzer = new QtlSeqAnalysis(this, options);
                return analyzer.Run();
            });
        }
    }
}
