namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// グラフ作成オプション
    /// </summary>
    public interface IGraphOptions
    {
        /// <summary>
        /// グラフ画像の幅(Pixel)を取得する。
        /// </summary>
        int FigureWidth { get; }

        /// <summary>
        /// グラフ画像の高さ(Pixel)を取得する。
        /// </summary>
        int FigureHeight { get; }

        /// <summary>
        /// X軸の目盛り間隔(MB)を取得する。
        /// </summary>
        int XAxisMajorStep { get; }
    }
}
